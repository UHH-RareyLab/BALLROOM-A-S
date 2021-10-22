import json
import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

import src.utils.definitions as definitions

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


class ScoreProcessor:
    def __init__(self):
        self.output_path = Path(definitions.output_path)

        # Initialize score data
        self.n_not_found = 0
        self.load_benchmark_files()
        self.load_scores()
        self.initialize_method_list()

        # Define cutoffs
        self.rmsd_cutoff = 2

        self.missing_score_dict = {x: set() for x in self.methods}

    def load_benchmark_files(self):
        """Load benchmark CSV files into pandas DataFrames."""
        self.benchmark_files = {
            "ballrooms": pd.read_csv(
                definitions.ballrooms / "complete_data.csv", index_col=0
            ),
            "ballrooma": pd.read_csv(
                definitions.ballrooma / "complete_data.csv", index_col=0
            ),
        }
        logging.info("Benchmarks loaded.")

    def load_scores(self):
        """Load all the score JSON files and CSVs into memory."""
        # Load JSON scores
        self.jamda_scores = self.load_json(definitions.score_path / "jamda_scores.json")
        self.boltz_affinity_scores = self.load_json(
            definitions.score_path / "boltz_pred_affinities.json"
        )
        self.boltz_pred_scores = self.load_json(
            definitions.score_path / "boltz_pred_probabilities.json"
        )
        # Load docking scores
        self.dock_scores = self.load_json(definitions.score_path / "DOCK.json")
        self.vina_scores = self.load_json(definitions.score_path / "Vina.json")
        # Load neural network output
        self.nn_scores_df = pd.read_csv(definitions.score_path / "nn_output.csv")
        self.nn_scores_dict = self.convert_scores_to_dict(
            self.nn_scores_df,
            key_target="target.chembl_id",
            key_ligand="Ligand",
            prediction="predictions",
        )
        self.jamda_rmsd = self.load_json(definitions.score_path / "jamda_rmsds.json")
        self.jamda_rescoring = self.load_json(
            definitions.score_path / "jamda_scores_rescoring.json"
        )
        logging.info("Scores loaded.")

    @staticmethod
    def load_json(path):
        """Load a JSON file with error handling."""
        try:
            with open(path, "r") as f:
                return json.load(f)
        except Exception as e:
            logging.error(f"Failed to load {path}: {e}")
            return {}

    @staticmethod
    def convert_scores_to_dict(df, key_target, key_ligand, prediction):
        """Convert DataFrame to nested dictionary for fast lookup."""
        result = {}
        for _, row in df.iterrows():
            target = row[key_target]
            ligand = row[key_ligand]
            pred = row[prediction]
            result.setdefault(target, {})[ligand] = pred
        return result

    def initialize_method_list(self):
        """Initialize the list of scoring methods."""
        self.methods = [
            "logp",
            "masse",
            "jamda",
            "tpsa",
            "boltz_pred",
            "boltz_affinity",
            "DOCK",
            "Vina",
            "nn_tobi",
            "jamda_rescoring",
            "jamda_workingrmsd",
        ]

    def get_score_from_dict(self, ligandname, dockinto, score_dict, method):
        """Retrieve the best score for a ligand from a score dictionary."""
        best_score = float("inf")
        for key in score_dict.keys():
            key = key.strip()
            for receptor in dockinto.split(";"):
                receptor = receptor.strip()
                if not receptor:
                    continue
                pdb = receptor.split("_")[0]
                ligand = receptor[5:]
                potential_key = f"{pdb}-{ligand}"
                if potential_key == key or pdb.lower() == key:
                    try:
                        if ligandname in score_dict[key]:
                            current_score = float(score_dict[key][ligandname])
                        else:
                            potential_keys = [
                                x
                                for x in score_dict[key].keys()
                                if x.startswith(ligandname)
                            ]
                            if not potential_keys:
                                raise KeyError("Ligand not found")
                            current_score = float("inf")
                            for pk in potential_keys:
                                score_value = float(score_dict[key][pk])
                                if score_value < current_score:
                                    current_score = score_value
                        if current_score < best_score:
                            best_score = current_score
                    except (KeyError, ValueError, TypeError) as e:
                        self.n_not_found += 1
                        logging.debug(
                            f"Score not found or invalid for method {method} and {ligandname}: {e}"
                        )
                        self.missing_score_dict[method].add(ligandname + receptor)
        if best_score > 1000:
            self.missing_score_dict[method].add(ligandname + receptor)
            logging.info(
                f"Score not found for method {method} and {ligandname} in {dockinto}, best_score={best_score}"
            )
        return best_score

    def get_both_scores(
        self, id1, id2, dockinto1, dockinto2, score_dict, factor=-1, method="test"
    ):
        """Get scores for both ligands, applying a factor."""
        score1 = factor * self.get_score_from_dict(
            id1, dockinto1, score_dict, method=method
        )
        score2 = factor * self.get_score_from_dict(
            id2, dockinto2, score_dict, method=method
        )
        return score1, score2

    def resolve_method_scores(self, score_metrics, method):
        if method in score_metrics:
            return score_metrics.get(method)
        if method.startswith("boltz_affinity"):
            return score_metrics.get("boltz_affinity")
        elif method.startswith("nn_tobi"):
            return score_metrics.get("nn_tobi")
        else:
            raise RuntimeError(f"Method {method} not known")

    def process_benchmark(self, name, df):
        """Process a specific benchmark DataFrame."""
        all_scores = {method: {"ges": [0, 0]} for method in self.methods}
        n_works = 0
        n_doesnt = 0

        for _, row in df.iterrows():
            # Extract data
            id1, id2 = row["ID1"], row["ID2"]
            smiles1, smiles2 = row["Smiles1"], row["Smiles2"]
            activity1, activity2 = float(row["Activity1"]), float(row["Activity2"])
            difference = int(row["Difference"])
            assay = row["Assay"]
            target = row["Target"]
            dockinto1, dockinto2 = row["Dockinto1"], row["Dockinto2"]

            # Convert SMILES to molecules
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if not mol1 or not mol2:
                logging.warning(f"Invalid SMILES for {id1} or {id2}")
                continue

            score_metrics = {}

            # Scores
            score_metrics["jamda"] = self.get_both_scores(
                id1, id2, dockinto1, dockinto2, self.jamda_scores, method="jamda"
            )
            score_metrics["boltz_affinity"] = self.get_both_scores(
                id1,
                id2,
                dockinto1,
                dockinto2,
                self.boltz_affinity_scores,
                method="boltz_affinity",
            )
            score_metrics["boltz_pred"] = self.get_both_scores(
                id1,
                id2,
                dockinto1,
                dockinto2,
                self.boltz_pred_scores,
                factor=1,
                method="boltz_pred",
            )
            score_metrics["Vina"] = self.get_both_scores(
                id1, id2, dockinto1, dockinto2, self.vina_scores, method="Vina"
            )
            score_metrics["DOCK"] = self.get_both_scores(
                id1, id2, dockinto1, dockinto2, self.dock_scores, method="DOCK"
            )

            # Additional descriptors
            score_metrics["logp"] = (Crippen.MolLogP(mol1), Crippen.MolLogP(mol2))
            score_metrics["masse"] = (
                Descriptors.ExactMolWt(mol1),
                Descriptors.ExactMolWt(mol2),
            )
            score_metrics["tpsa"] = (Descriptors.TPSA(mol1), Descriptors.TPSA(mol2))

            # NN scores
            score_metrics["nn_tobi"] = (
                self.nn_scores_dict.get(target, {}).get(id1, None),
                self.nn_scores_dict.get(target, {}).get(id2, None),
            )
            score_metrics["jamda_rescoring"] = self.get_both_scores(
                id1,
                id2,
                dockinto1,
                dockinto2,
                self.jamda_rescoring,
                method="jamda_rescoring",
            )
            score_metrics["jamda_workingrmsd"] = self.get_both_scores(
                id1,
                id2,
                dockinto1,
                dockinto2,
                self.jamda_scores,
                method="jamda_workingrmsd",
            )

            # Compare scores and update metrics
            for method in self.methods:
                val1 = self.resolve_method_scores(score_metrics, method)
                # val1 = score_metrics.get(method)
                if val1 is None:
                    continue
                score1, score2 = val1

                # Apply filtering based on cutoffs
                if "cliffs" in method and difference < 2:
                    continue
                if method == "jamda_workingrmsd":
                    if id1 not in self.jamda_rmsd or id2 not in self.jamda_rmsd:
                        # Not all rmsds present
                        print(
                            f"RMSD not found for {id1} or {id2} {id1 in self.jamda_rmsd}, {id2 in self.jamda_rmsd}"
                        )
                        continue
                    else:
                        if (
                            self.jamda_rmsd[id1] > self.rmsd_cutoff
                            or self.jamda_rmsd[id2] > self.rmsd_cutoff
                        ):
                            print(
                                f"Threshhold not met for {id1}{self.jamda_rmsd[id1]} or {id2}{self.jamda_rmsd[id2]}"
                            )
                            # One or more RMSDs bad
                            continue

                # Decide which score indicates higher activity
                if activity1 < activity2:
                    if score1 > score2:
                        self.increment_score(all_scores, method, target)
                        n_works += 1
                    else:
                        self.increment_score(all_scores, method, target, success=False)
                        n_doesnt += 1
                else:
                    if score2 > score1:
                        self.increment_score(all_scores, method, target)
                        n_works += 1
                    else:
                        self.increment_score(all_scores, method, target, success=False)
                        n_doesnt += 1

        # Save results
        self.save_results(name, all_scores)
        return n_works, n_doesnt

    def increment_score(self, all_scores, method, target, success=True):
        """Increment positive/total counts."""
        if target not in all_scores[method]:
            all_scores[method][target] = [0, 0]
        if success:
            all_scores[method][target][0] += 1
            all_scores[method]["ges"][0] += 1
        all_scores[method][target][1] += 1
        all_scores[method]["ges"][1] += 1

    def save_results(self, benchmark_name, data):
        """Save the results to CSV and JSON."""
        # Prepare DataFrame
        records = []
        for method, targets in data.items():
            for target, counts in targets.items():
                positives, total = counts
                ratio = round(positives / total, 4) * 100 if total else 0
                records.append([method, target, positives, total, ratio])
        df = pd.DataFrame(
            records, columns=["Method", "Target", "Positives", "Total", "Ratio"]
        )
        df.to_csv(self.output_path / f"{benchmark_name}_rmsd.csv", index=False)
        with open(self.output_path / f"{benchmark_name}_rmsd.json", "w") as f:
            json.dump(data, f)

    def run(self):
        """Main execution: process all benchmarks."""
        for name, df in self.benchmark_files.items():
            logging.info(f"Processing {name}...")
            n_works, n_doesnt = self.process_benchmark(name, df)
            logging.info(f"{name} completed. Success: {n_works}, Failures: {n_doesnt}")
            logging.info(f"Not found score count: {self.n_not_found}")


if __name__ == "__main__":

    processor = ScoreProcessor()
    processor.run()
