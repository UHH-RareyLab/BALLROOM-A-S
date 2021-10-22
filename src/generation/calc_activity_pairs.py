import json
import shutil

import pandas as pd

from src.generation.check_pair_benchmark_contradictions import get_noncontradicting_set
from src.utils.definitions import pair_data_columns


def calc_difference_between_two_activities(activity1: float, activity2: float):
    activity1 = float(activity1)
    activity2 = float(activity2)
    i = 0
    found = False
    while not found:
        i += 1
        if activity1 > activity2:
            if pow(10, i) * activity2 > activity1:
                break
        else:
            if pow(10, i) * activity1 > activity2:
                break
        if i > 100:
            raise RuntimeError(i, activity1, activity2)
    return i - 1


class ActivityPair:
    def __init__(
        self,
        identifier1: str,
        identifier2: str,
        smiles1: str,
        smiles2: str,
        activity1: float,
        activity2: float,
        difference: int,
        assay: str,
        target: str,
        activity_type: str,
        activity_unit: str,
        dock_into1: str,
        dock_into2: str,
        activity_id_1: str,
        activity_id_2: str,
    ):
        self.id1 = identifier1
        self.id2 = identifier2
        self.smiles1 = smiles1
        self.smiles2 = smiles2
        self.activity1 = activity1
        self.activity2 = activity2
        self.assay = assay
        self.target = target
        self.difference = difference
        self.activity_type = activity_type
        self.activity_unit = activity_unit
        self.dock_into1 = dock_into1
        self.dock_into2 = dock_into2
        self.activity_id_1 = activity_id_1
        self.activity_id_2 = activity_id_2

    def report_pair(self):
        print(f"IDs: {self.id1} {self.id2} Smiles: {self.smiles1} {self.smiles2}")
        print(f"Target {self.target} Assay {self.assay}")
        print(
            f"Activities {self.activity1} {self.activity2} {self.difference} {self.dock_into1} {self.activity_id_1} {self.activity_id_2}"
        )

    def return_pair_as_list(self):
        return [
            self.id1,
            self.id2,
            self.smiles1,
            self.smiles2,
            self.activity1,
            self.activity2,
            self.difference,
            self.activity_type,
            self.activity_unit,
            self.assay,
            self.target,
            self.dock_into1,
            self.dock_into2,
            self.activity_id_1,
            self.activity_id_2,
        ]


def calc_pairs_mmp_assay(assay_path, target_name, assay_name, min_difference=1):
    new_pairs = []
    assay = assay_name
    target = target_name
    all_pdb_names = ""
    data = pd.read_csv(assay_path)
    for index, row in data.iterrows():
        identifier1 = row["Name"]
        if not identifier1.startswith("CHEMBL"):
            all_pdb_names = all_pdb_names + identifier1 + ";"
    if all_pdb_names == "":
        return new_pairs
    for index1, row1 in data.iterrows():
        for index2, row2 in data.iterrows():
            if index1 >= index2:
                continue
            affinity1 = row1["Activity"]
            affinity2 = row2["Activity"]
            identifier1 = row1["Name"]
            identifier2 = row2["Name"]
            smiles1 = row1["SMILES"]
            smiles2 = row2["SMILES"]
            activity_type = row1["Activity Type"]
            activity_unit = row1["Activity Unit"]
            activity_id1 = row1["Activity_id"]
            activity_id2 = row2["Activity_id"]
            difference = calc_difference_between_two_activities(affinity1, affinity2)
            if difference >= min_difference:
                new_pairs.append(
                    ActivityPair(
                        identifier1,
                        identifier2,
                        smiles1,
                        smiles2,
                        affinity1,
                        affinity2,
                        difference,
                        assay,
                        target,
                        activity_type,
                        activity_unit,
                        all_pdb_names,
                        all_pdb_names,
                        activity_id1,
                        activity_id2,
                    )
                )
    return new_pairs


def calc_pairs_in_assay(assay_path, min_difference=1):
    new_pairs = []
    data = pd.read_csv(assay_path)
    for index1, row1 in data.iterrows():
        assay = row1["assay.chembl_id"]
        target = row1["target.chembl_id"]
        for index2, row2 in data.iterrows():
            if index1 >= index2:
                continue
            affinity1 = row1["activity.standard_value"]
            affinity2 = row2["activity.standard_value"]
            identifier1 = row1["pdb.pdb_id"] + "_" + row1["name"]
            identifier2 = row2["pdb.pdb_id"] + "_" + row2["name"]
            smiles1 = row1["usmiles"]
            smiles2 = row2["usmiles"]
            activity_type = row1["activity.standard_type"]
            activity_unit = row1["activity.standard_units"]
            activity_id1 = row1["activity.chembl_activity_id"]
            activity_id2 = row2["activity.chembl_activity_id"]
            difference = calc_difference_between_two_activities(affinity1, affinity2)
            if difference >= min_difference:
                new_pairs.append(
                    ActivityPair(
                        identifier1,
                        identifier2,
                        smiles1,
                        smiles2,
                        affinity1,
                        affinity2,
                        difference,
                        assay,
                        target,
                        activity_type,
                        activity_unit,
                        identifier1,
                        identifier2,
                        activity_id1,
                        activity_id2,
                    )
                )
    return new_pairs


def create_pair_benchmark_crystal(input_path, output_path, min_difference):
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir()
    new_benchmark_data = {}
    for target in input_path.glob("*"):
        if target.is_file():
            continue
        new_benchmark_data[target.name] = []
        for assay in target.glob("*.csv"):
            assay_pairs = calc_pairs_in_assay(assay, min_difference)
            new_benchmark_data[target.name].extend(assay_pairs)
    complete_pairs = 0
    n_targets = 0
    complete_data = []
    contradiction_data = []
    for target in new_benchmark_data:
        if len(new_benchmark_data[target]) == 0:
            continue
        current_target_path = output_path / target
        current_target_path.mkdir()
        for assay_data in (input_path / target).glob("*.csv"):
            shutil.copyfile(assay_data, current_target_path / assay_data.name)
        activity_pair_data = []
        for pair in new_benchmark_data[target]:
            pair_data = pair.return_pair_as_list()
            activity_pair_data.append(pair_data)
        df = pd.DataFrame(
            activity_pair_data,
            columns=pair_data_columns,
        )
        (
            contradicting_assays,
            _,
            _,
            _,
            assay_labels,
            contradiction_matrix,
        ) = get_noncontradicting_set(df, target)
        print(
            f"contradicting assays for target {target} in crystal benchmark:",
            contradicting_assays,
        )
        contradiction_data.append([target, len(contradicting_assays)])
        activity_pair_data_true = [
            x for x in activity_pair_data if x[9] not in contradicting_assays
        ]
        complete_pairs += len(activity_pair_data_true)
        df_true = pd.DataFrame(
            activity_pair_data_true,
            columns=pair_data_columns,
        )
        if len(contradicting_assays) > 0:
            with open(current_target_path / "contradicting_assays.json", "w") as f:
                json.dump(contradicting_assays, f)
            contradiction_matrix = pd.DataFrame(
                contradiction_matrix, index=assay_labels, columns=assay_labels
            )
            contradiction_matrix.to_csv(
                current_target_path / "contradiction_matrix.csv"
            )
        if len(activity_pair_data_true) > 0:
            complete_data.extend(activity_pair_data_true)
            n_targets += 1
            df_true.to_csv(current_target_path / "pair_data.csv")
        else:
            shutil.rmtree(current_target_path)

    complete_data = pd.DataFrame(complete_data, columns=pair_data_columns)
    complete_data.to_csv(output_path / "complete_data.csv")
    contradiction_data = pd.DataFrame(
        contradiction_data, columns=["Target", "N Contradictions"]
    )
    contradiction_data.to_csv(output_path / "contradiction_data.csv")


def create_pair_benchmark_mmp(input_path, output_path, min_difference):
    if output_path.exists():
        shutil.rmtree(output_path)
    output_path.mkdir()
    new_benchmark_data = {}
    for target in input_path.glob("*"):
        if target.is_file():
            continue
        new_benchmark_data[target.name] = []
        for assay in target.glob("*"):
            if assay.is_file():
                continue
            for mmp in assay.glob("*.csv"):
                assay_pairs = calc_pairs_mmp_assay(
                    mmp, target.name, assay.name.split("-")[0], min_difference
                )
                current_target_path = output_path / target.name
                current_assay_path = output_path / target.name / assay.name
                if len(assay_pairs) > 0:
                    if not current_target_path.exists():
                        current_target_path.mkdir()
                    if not current_assay_path.exists():
                        current_assay_path.mkdir()
                    shutil.copyfile(mmp, current_assay_path / mmp.name)
                    shutil.copyfile(
                        mmp.parent / (mmp.stem + "_data.json"),
                        current_assay_path / (mmp.stem + "_data.json"),
                    )
                    shutil.copyfile(
                        mmp.parent / (mmp.stem + ".smi"),
                        current_assay_path / (mmp.stem + ".smi"),
                    )
                    new_benchmark_data[target.name].extend(assay_pairs)
    complete_pairs = 0
    n_targets = 0
    complete_data = []
    contradiction_data = []
    for target in new_benchmark_data:
        if len(new_benchmark_data[target]) == 0:
            continue
        current_target_path = output_path / target
        activity_pair_data = []
        for pair in new_benchmark_data[target]:
            pair_data = pair.return_pair_as_list()
            activity_pair_data.append(pair_data)
        df = pd.DataFrame(
            activity_pair_data,
            columns=pair_data_columns,
        )
        (
            contradicting_assays,
            _,
            _,
            _,
            assay_labels,
            contradiction_matrix,
        ) = get_noncontradicting_set(df, target)
        print(
            f"contradicting assays for target {target} in mmp benchmark:",
            contradicting_assays,
        )
        activity_pair_data_true = [
            x for x in activity_pair_data if x[9] not in contradicting_assays
        ]
        df_true = pd.DataFrame(
            activity_pair_data_true,
            columns=pair_data_columns,
        )
        contradiction_data.append([target, len(contradicting_assays)])
        if len(contradicting_assays) > 0:
            with open(current_target_path / "contradicting_assays.json", "w") as f:
                json.dump(contradicting_assays, f)
            contradiction_matrix = pd.DataFrame(
                contradiction_matrix, index=assay_labels, columns=assay_labels
            )
            contradiction_matrix.to_csv(
                current_target_path / "contradiction_matrix.csv"
            )
        if len(activity_pair_data_true) > 0:
            complete_data.extend(activity_pair_data_true)
            n_targets += 1
            complete_pairs += len(activity_pair_data_true)
            df_true.to_csv(current_target_path / "pair_data.csv")
        else:
            shutil.rmtree(current_target_path)
    complete_data = pd.DataFrame(complete_data, columns=pair_data_columns)
    complete_data.to_csv(output_path / "complete_data.csv")
    contradiction_data = pd.DataFrame(
        contradiction_data, columns=["Target", "N Contradictions"]
    )
    contradiction_data.to_csv(output_path / "contradiction_data.csv")
