import json
import math
import shutil
from pathlib import Path
from typing import Dict, List

import pandas as pd
import progressbar

from src.generation.datasetFilter import (
    DatasetFilter,
    allowed_activities,
    matchings_proteins,
    matchings_molecules,
)
from src.generation.download_benchmark import download_benchmark
from src.generation.finalize_siena_benchmark import (
    finalize_siena_benchmark,
)
from src.generation.make_siena_clusters import make_siena_clusters
from src.generation.retrieve_bioactivity_info_from_chembl import (
    download_additional_assay_data_for_benchmark
)
from src.utils.definitions import (
    mm_strict,
    siena_path,
    sienatools_path,
    stractablef_path,
    general_filters,
)


class ActivityBenchmarkMaker:
    def run(
        self,
        min_datapoints: int,
        min_difference: int,
        ediam: float,
        benchmark_name: str = "affinity_benchmark",
        allowed_protein_matches: List[str] = matchings_proteins,
        allowed_molecule_matches: List[str] = matchings_molecules,
        only_preliminary: bool = False,
        filter_uncertain=True,
        max_mutations=0,
        pdbs_included_in_assays=[],
        pdbs_included_total=[],
        targets_allowed=[],
    ):
        """
        This function initiates the affinity benchmark creation
        After calculating all final dataframes, it exports the dataframes to harddrive
        @param min_datapoints: Minimum number of datapoints present in a testcase for it to be added to the benchmark
        @param min_difference: Mimimum difference in oder of magnitude in a single testcase
        @param ediam: Minumum EDIAm allowed for molecules
        @param benchmark_name: Name of the folder in which the benchmark is exported
        @param allowed_protein_matches : List of allowed protein matchings
        @param allowed_molecule_matches : List of allowed molecule matchings
        @param pdbs_included_in_assays : List of pdbs that an assay has to have at least one of to be exported
        @param pdbs_included_total: List of pdbs that are not discarded upen first readin
        @return: None
        """
        basic_filtered_data = []
        i = 0
        n_pass = 0
        n_notpass = 0
        print("Reading in the master table")
        with open(stractablef_path, "r") as f:
            length = sum(1 for _ in f)
        chunksize = 10000
        steps = int(length / chunksize)
        widget = [
            "Loading Master Table: ",
            progressbar.BouncingBar(),
            progressbar.AdaptiveETA(),
        ]
        for chunk in progressbar.progressbar(
            pd.read_csv(
                stractablef_path,
                index_col=0,
                chunksize=chunksize,
                dtype=str,
                skiprows=1,
            ),
            max_value=steps,
            widgets=widget,
        ):
            columns = chunk.columns.values
            for index, row in chunk.iterrows():
                i += 1
                worked = self.basic_filter_row(
                    row,
                    allowed_protein_matches,
                    allowed_molecule_matches,
                    ediam,
                    max_mutations,
                    pdbs_included_total,
                    targets_allowed,
                )
                if worked:
                    basic_filtered_data.append([x for x in row])
                    n_pass += 1
                else:
                    n_notpass += 1
        print("Filtering the data")
        basic_filtered_data = pd.DataFrame(
            basic_filtered_data, columns=columns, dtype=str
        )
        print(i, basic_filtered_data.shape)
        final_targets = self.calc_targets_assays(basic_filtered_data, min_datapoints)
        initial_dataframes = self.make_preliminary_dataset(
            final_targets, basic_filtered_data
        )
        print("Create preliminary benchmark")
        preliminary_dataframes = self.filter_dataframes(
            preliminary_dataframes=initial_dataframes,
            min_datapoints=min_datapoints,
            min_difference=min_difference,
            full_filtering=False,
            filter_uncertain=filter_uncertain,
            pdbs_included_in_assay=pdbs_included_in_assays,
        )
        print("Export preliminary dataset")
        big_data = self.export_dataset(
            preliminary_dataframes, benchmark_name + "_preliminary"
        )
        self.export_settings_results(
            min_datapoints,
            min_difference,
            ediam,
            allowed_protein_matches,
            allowed_molecule_matches,
            big_data[0],
            big_data[1],
            big_data[2],
            benchmark_name + "_preliminary",
        )
        download_benchmark(Path().cwd() / (benchmark_name + "_preliminary"))
        if only_preliminary:
            return
        make_siena_clusters(
            path_benchmark=Path().cwd() / (benchmark_name + "_preliminary"),
            path_siena=siena_path,
            path_sienatools=sienatools_path,
        )
        print("Finalize Siena benchmark")
        all_dataframes = finalize_siena_benchmark(
            benchmark_path=Path().cwd() / (benchmark_name + "_preliminary"),
            min_datapoints=min_datapoints,
        )
        print("Final Checkup and Export")
        final_dataframes = self.filter_dataframes(
            preliminary_dataframes=all_dataframes,
            min_datapoints=min_datapoints,
            min_difference=min_difference,
            full_filtering=True,
        )
        small_data = self.export_dataset(
            final_dataframes=final_dataframes,
            foldername=benchmark_name,
            triplename=True,
        )
        self.export_settings_results(
            min_datapoints,
            min_difference,
            ediam,
            allowed_protein_matches,
            allowed_molecule_matches,
            small_data[0],
            small_data[1],
            small_data[2],
            benchmark_name,
        )
        benchmark_path = Path().cwd() / benchmark_name
        download_additional_assay_data_for_benchmark(benchmark_path)

    def export_settings_results(
        self,
        min_datapoints,
        min_difference,
        ediam,
        protein_matchings,
        ligand_matchings,
        number_targets,
        number_assays,
        number_rows,
        foldername,
    ):
        settings_results = {
            "Foldername": foldername,
            "min_datapoints": min_datapoints,
            "min_difference": min_difference,
            "protein_matchings": protein_matchings,
            "ligand_matchings": ligand_matchings,
            "ediam": ediam,
            "number_targets": number_targets,
            "number_assays": number_assays,
            "number_rows": number_rows,
        }
        jsonstring = json.dumps(settings_results, indent=4)
        print(jsonstring)
        with open(Path(foldername) / "settings_results.json", "w") as outfile:
            outfile.write(jsonstring)

    def basic_filter_row(
        self,
        row,
        allowed_protein_matches,
        allowed_ligand_matches,
        ediam,
        max_mutations,
        pdbs_included,
        allowed_targets,
    ):
        """
        This function checks if a row of the master table adheres to the filter criteria
        @param row: row to check
        @param allowed_protein_matches: allowed matching levels of the protein matching
        @param allowed_ligand_matches: allowed matching levels of the molecule matching
        @param ediam: minumum ediam allowed
        @return: True if row passed, false if not
        """
        if (
            row["first_merge_ligandextractor_structureprofiler"] != "both"
            or row["second_merge_activityfinder_rest"] != "both"
            or row["skip_reason"] != "NotSkipped"
            or row["activity.standard_relation"] != "="
            or float(row["EDIAm"]) < ediam
            or math.isnan(float(row["EDIAm"]))
            or not row["activity.standard_type"] in allowed_activities
            or int(row["heavy_atoms"]) < general_filters["heavy_atoms_min"]
            or int(row["heavy_atoms"]) > general_filters["heavy_atoms_max"]
            or not row["Level_pm"] in allowed_protein_matches
            or not row["small_molecule_matching_confidence_level.comment"]
            in allowed_ligand_matches
            or float(row["resolution"]) > general_filters["resolution_max"]
            or float(row["DPI"]) > general_filters["dpi_max"]
            or float(row["rFactor"]) > general_filters["rfactor_max"]
            or row["overfittingTest"] != "True"
            or row["significanceTest"] != "True"
            or row["noCrystalContacts"] != "True"
            or int(row["N_mutations"]) > max_mutations
            or (
                row["target.chembl_id"] not in allowed_targets
                and len(allowed_targets) > 0
            )
            or (row["pdb"] not in pdbs_included and len(pdbs_included) > 0)
        ):
            return False
        else:
            return True

    @staticmethod
    def filter_dataframes(
        preliminary_dataframes: Dict[tuple, pd.DataFrame],
        min_datapoints,
        min_difference,
        full_filtering=False,
        filter_uncertain=True,
        pdbs_included_in_assay=[],
    ) -> Dict[tuple, pd.DataFrame]:
        """
        Filters the testcases so that they conform to the standards of this benchmark
        @param preliminary_dataframes: The dataframes of testcases to be filtered
        @param min_datapoints: Minimum number of datapoints present in a testcase for it to be added to the benchmark
        @param min_difference: Mimimum difference in oder of magnitude in a single testcase
        @return: The filtered dataframes
        """
        final_dataframes = {}
        Filter = DatasetFilter()
        for dataset in preliminary_dataframes:
            if len(pdbs_included_in_assay) > 0:
                pdbs_in_dataset = [
                    x.lower() for x in preliminary_dataframes[dataset]["pdb"]
                ]
                if (
                    len([x for x in pdbs_in_dataset if x in pdbs_included_in_assay])
                    == 0
                ):
                    continue
            try:
                if full_filtering:
                    final_dataframes[dataset] = Filter.run(
                        preliminary_dataframes[dataset], min_datapoints, min_difference
                    )
                else:
                    final_dataframes[dataset] = Filter.run_prefiltering(
                        preliminary_dataframes[dataset],
                        min_datapoints,
                        filter_uncertain=filter_uncertain,
                    )
            except RuntimeError:
                continue
            # except ValueError:
            #    print("ERROR HERE LOOK LISTEN")
            #    print(dataset)
            #    continue
        return final_dataframes

    def make_smi_from_pandas(self, dataframe: pd.DataFrame, full_path: Path):
        """
        Makes a smiles file from a pandas dataframe (basically just a tsv)
        @param dataframe: The dataframe that it shall use
        @param full_path: The path to which it shall save the smiles file
        """
        molecule_data = []
        for index, row in dataframe.iterrows():
            usmiles = row["usmiles"]
            pdb = row["pdb"]
            name = row["name"]
            molecule_data.append([usmiles, pdb + "_" + name])
        molecule_data = pd.DataFrame(molecule_data, dtype=str)
        molecule_data.to_csv(full_path, sep=" ", header=False, index=False)

    def export_dataset(
        self,
        final_dataframes: Dict[tuple, pd.DataFrame],
        foldername: str,
        triplename=False,
    ):
        """
        Exports all dataframes in a dict
        @param final_dataframes: The dict containing the dataframes to be exported
        @param foldername: the name of the folder that the benchmark shall be exported in
        @param triplename: State if the tuple of the dict is consisting of two or three parts
        @return: None
        """
        full_dataframe = pd.DataFrame()
        foldername = Path(foldername)
        if foldername.is_dir():
            shutil.rmtree(foldername)
        foldername.mkdir()
        n_assays = 0
        n_targets = 0
        n_rows = 0
        for dataset in final_dataframes:
            to_export = final_dataframes[dataset]
            base_path = foldername / dataset[0]
            if not base_path.is_dir():
                n_targets += 1
                base_path.mkdir()
            if triplename:
                savename = base_path / f"{dataset[1]}-{dataset[2]}"
            else:
                savename = base_path / f"{dataset[1]}"
            if (Path(str(savename) + ".csv")).exists():
                continue
            to_export.to_csv(str(savename) + ".csv")
            full_dataframe = pd.concat([full_dataframe, to_export])
            self.make_smi_from_pandas(to_export, Path(str(savename) + ".smi"))
            n_assays += 1
            n_rows += to_export.shape[0]
        print(
            f"Number of targets: {n_targets} and Number of assays correctly exported: {n_assays} with a total of {n_rows} rows"
        )
        full_dataframe.to_csv(foldername / "complete_benchmark.csv")
        return [n_targets, n_assays, n_rows]

    def make_preliminary_dataset(
        self,
        final_target_map: Dict[tuple, set],
        dataframe: pd.DataFrame,
        pure_activities=False,
    ) -> Dict[tuple, pd.DataFrame]:
        """
        Creates the preliminary dataset using
        @param final_target_map: Dict using target/assay tuple and ligand Usmiles as set
        @param dataframe: complete dataframe that shall be filtered
        @return: Dict using target/assay tuple and dataframe part of the whole dataframe as value
        """
        final_datasets = {}
        for index, row in dataframe.iterrows():
            target = row["target.chembl_id"]
            if pure_activities:
                ligand = row["pdb_ligand.naomi_unique_smiles"]
            else:
                ligand = row["usmiles"]
            assay = row["assay.chembl_id"]
            for key in final_target_map.keys():
                if (
                    target == key[0]
                    and assay == key[1]
                    and ligand in final_target_map[key]
                ):
                    if key in final_datasets.keys():
                        final_datasets[key].append(row)
                    else:
                        final_datasets[key] = [row]

        for dataset in final_datasets:
            final_datasets[dataset] = pd.DataFrame(final_datasets[dataset], dtype=str)
        return final_datasets

    def calc_targets_assays(
        self, dataframe: pd.DataFrame, min_observations, pure_activities=False
    ) -> Dict[tuple, set]:
        """
        Calculates all target assay combinations present in the dataframe
        Checks that there are at least min_observations ligands in the dataframe to not overconvolute the target/assays
        @param dataframe: Raw bob dataframe
        @return: Dict with tuple target/assay as key and set of ligand smiles as value
        """
        map_target_ligand_assay = {}

        for index, row in dataframe.iterrows():
            target = row["target.chembl_id"]
            if pure_activities:
                ligand = row["pdb_ligand.naomi_unique_smiles"]
            else:
                ligand = row["usmiles"]
            assays = row["assay.chembl_id"]
            for assay in (
                assays.replace("[", "").replace("]", "").replace("'", "").split(",")
            ):
                if target in map_target_ligand_assay.keys():
                    if assay in map_target_ligand_assay[target].keys():
                        map_target_ligand_assay[target][assay].add((ligand))
                    else:
                        map_target_ligand_assay[target][assay] = {ligand}
                else:
                    map_target_ligand_assay[target] = {assay: {ligand}}

        final_targets_assays = {}
        for target in map_target_ligand_assay.keys():
            for assay in map_target_ligand_assay[target].keys():
                n_ligands_assay = len(map_target_ligand_assay[target][assay])
                all_ligands = map_target_ligand_assay[target][assay]
                if n_ligands_assay >= min_observations:
                    key = (target, assay)
                    final_targets_assays[key] = all_ligands
        return final_targets_assays

    def custom_complex_filter(self, dataset: pd.DataFrame) -> pd.DataFrame:
        """
        This is the custom complex filter that is used in the WIP redocking the PDB paper and was taken from
        the StructureProfiler Paper
        @param dataset: dataset to be filtered
        @return: filtered dataset
        """
        new_dataset = dataset[dataset["resolution"] <= 2.5]
        new_dataset = new_dataset[new_dataset["DPI"] <= 0.42]
        new_dataset = new_dataset[new_dataset["rFactor"] <= 0.42]
        new_dataset = new_dataset[new_dataset["rFree"] <= 0.45]
        new_dataset = new_dataset[new_dataset["overfittingTest"] == True]
        new_dataset = new_dataset[new_dataset["significanceTest"] == True]
        return new_dataset


def main():
    benchmark = ActivityBenchmarkMaker()
    benchmark.run(
        min_datapoints=2,
        min_difference=1,
        ediam=0.4,
        allowed_protein_matches=["Platin", "Gold", "Silber", "Bronze"],
        allowed_molecule_matches=mm_strict,
        benchmark_name="pair_benchmark_starting_point_strict_mols",
    )


if __name__ == "__main__":
    main()
