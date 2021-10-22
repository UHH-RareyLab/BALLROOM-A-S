import shutil
from itertools import product
from pathlib import Path

import Levenshtein
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components

from src.utils import db_interaction as db


def find_bioactivities_for_assay_alt(assay):
    query = (
        f"SELECT a.activity_id, a.data_validity_comment, a.potential_duplicate, a.standard_relation, a.standard_value, a.standard_type, m.chembl_id as molecule_chembl_id, c.canonical_smiles, t.chembl_id as target_chembl_id FROM assays "
        f"JOIN public.activities a on assays.assay_id = a.assay_id "
        f"JOIN public.molecule_dictionary m on a.molregno = m.molregno "
        f"JOIN public.compound_structures c on a.molregno = c.molregno "
        f"JOIN target_dictionary t on assays.tid = t.tid "
        f"WHERE assays.chembl_id = '{assay}'"
    )
    return db.query_chembl(query)


# ====================================================================#
def check_validity_comment(comment):
    allowed = [
        "Manually validated",
        "Potential missing data",
        "Potential author error",
        "Potential transcription error",
    ]
    if pd.isna(comment) or comment in allowed:
        return True
    return False


def get_activities_smiles(assay_id, activity_type):
    raw_data = find_bioactivities_for_assay_alt(assay_id)
    final_results = []
    for index, activity in raw_data.iterrows():
        data_validity_comment = activity["data_validity_comment"]
        duplicate = activity["potential_duplicate"]
        relation = activity["standard_relation"]
        molecule_chembl_id = activity["molecule_chembl_id"]
        type = activity["standard_type"]
        value = activity["standard_value"]
        smiles = activity["canonical_smiles"]
        activity_id = activity["activity_id"]
        if (
            check_validity_comment(data_validity_comment)
            and duplicate == 0
            and type == activity_type
            and relation == "="
        ):
            final_results.append([smiles, value, molecule_chembl_id, activity_id])
    return final_results


def get_data_for_assay(assay_id):
    raw_data = find_bioactivities_for_assay_alt(assay_id)
    activity_dict = {}
    for index, activity in raw_data.iterrows():
        activity_id = activity["activity_id"]
        target_id = activity["target_chembl_id"]
        molecule_chembl_id = activity["molecule_chembl_id"]
        data_validity_comment = activity["data_validity_comment"]
        duplicate = activity["potential_duplicate"]
        activity_dict[activity_id] = [
            check_validity_comment(data_validity_comment) and duplicate == 0,
            check_validity_comment(data_validity_comment),
            duplicate == 0,
            data_validity_comment,
        ]
    return activity_dict


def get_activities_for_assay(assay_id):
    raw_data = find_bioactivities_for_assay_alt(assay_id)
    activity_dict = {}
    for index, activity in raw_data.iterrows():
        activity_id = activity["activity_id"]
        target_id = activity["target_chembl_id"]
        molecule_chembl_id = activity["molecule_chembl_id"]
        data_validity_comment = activity["data_validity_comment"]
        duplicate = activity["potential_duplicate"]
        activity_dict[activity_id] = (
            check_validity_comment(data_validity_comment) and duplicate == 0
        )
    return activity_dict


def get_components_for_target(target_id):
    query = (
        f"SELECT c.component_id from target_dictionary t "
        f"JOIN target_components c on t.tid = c.tid "
        f"WHERE t.chembl_id = '{target_id}'"
    )
    data = db.query_chembl(query)
    return data["component_id"].tolist()


def get_sequence_for_component_id(component_id):
    query = f"SELECT sequence from component_sequences WHERE component_sequences.component_id = '{component_id}'"
    data = db.query_chembl(query)
    return data["sequence"].tolist()


def get_target_sequences_dict(target_list):
    target_sequences = {}
    for target in target_list:
        data = get_components_for_target(target)
        if len(data) == 0:
            target_sequences[target] = []
            continue
        for component in data:
            sequence = get_sequence_for_component_id(component)
            if target in target_sequences.keys():
                target_sequences[target].append(sequence)
            else:
                target_sequences[target] = [sequence]
    return target_sequences


def get_target_clusters(target_list, cutoff=0.95):
    target_sequences = get_target_sequences_dict(target_list)
    clustermatrix = create_clusters(target_sequences, cutoff)
    number_components, order_components = connected_components(clustermatrix)
    merged_target_list = [[] for _ in range(number_components)]
    for index, target in enumerate(target_list):
        merged_target_list[order_components[index]].append(target)
    return merged_target_list


def calc_levenshtein_two_targets(sequences1, sequences2):
    best_match = 0
    output = [x for x in product(sequences1, sequences2)]
    for pair in output:
        seq_identity = Levenshtein.ratio(pair[0], pair[1])
        if seq_identity > best_match:
            best_match = seq_identity
    return best_match


def initiate_adjacency_matrix(n):
    """
    Initiates an symetric adjacency matrix of size n with diagonal at 1 and rest at 0
    @param n: size of adjacency matrix
    @return: initiated adjacency matrix
    """
    matrix = np.full((n, n), 0)
    for i in range(n):
        matrix[i][i] = 1
    return matrix


def create_clusters(target_sequences, cutoff):
    adjacency_matrix = initiate_adjacency_matrix(len(target_sequences.keys()))
    for index1, target in enumerate(target_sequences):
        for index2, target2 in enumerate(target_sequences):
            if index1 <= index2:
                continue
            sequences1 = target_sequences[target]
            sequences2 = target_sequences[target2]
            best_match = calc_levenshtein_two_targets(sequences1, sequences2)
            works = best_match >= cutoff
            adjacency_matrix[index1][index2] = works
            adjacency_matrix[index2][index1] = works
    return adjacency_matrix


def merge_and_copy(merged_target_list, benchmark_path):
    new_benchmark = benchmark_path.parent / (benchmark_path.name + "_merged")
    if new_benchmark.exists():
        shutil.rmtree(new_benchmark)
    new_benchmark.mkdir()
    for files in benchmark_path.glob("*"):
        if files.is_file():
            shutil.copyfile(files, new_benchmark / files.name)
    for cluster in merged_target_list:
        if len(cluster) == 1:
            shutil.copytree(benchmark_path / cluster[0], new_benchmark / cluster[0])
        else:
            new_name = "_".join(cluster)
            new_folder = new_benchmark / new_name
            new_folder.mkdir()
            for target in cluster:
                for file in (benchmark_path / target).glob("*"):
                    shutil.copyfile(file, (new_folder / file.name))


def download_additional_assay_data_for_benchmark(benchmark_path):
    for target in benchmark_path.glob("*"):
        if target.is_file():
            continue
        for assay in target.glob("*.csv"):
            existing_mols = []
            data = pd.read_csv(assay, index_col=0)
            assay_type = data.iloc[0]["activity.standard_type"]
            for index, row in data.iterrows():
                existing_mols.append(row["chembl_ligand.chembl_id"])
            maybe_new_molecules = get_activities_smiles(assay.stem[:-2], assay_type)
            actually_new_molecules = [
                x
                for x in maybe_new_molecules
                if (not x[2] in existing_mols) and (not x[0] == None)
            ]
            new_mols_for_smiles = [[x[0], x[2]] for x in actually_new_molecules]
            if len(actually_new_molecules) == 0:
                continue
            actually_new_molecules = pd.DataFrame(
                actually_new_molecules,
                columns=["Smiles", "Activity", "ChemblID", "Activity_id"],
            )
            new_mols_for_smiles = pd.DataFrame(new_mols_for_smiles)
            actually_new_molecules.to_csv(
                assay.parent / (assay.stem + "_additional.tsv"), sep="\t"
            )
            new_mols_for_smiles.to_csv(
                assay.parent / (assay.stem + "_additional.smiles"),
                sep="\t",
                header=False,
                index=False,
            )


def main():
    data = find_bioactivities_for_assay_alt("CHEMBL839790")
    new_data = get_data_for_assay("CHEMBL839790")
    assay_path = Path(
        "/work/gutermuth/bob_analyse_skripte/gutermuth/benchmark_directory/presubmission_activityfinder_bugfix/raw_benchmark_start_merged/CHEMBL239/CHEMBL839790-0.csv"
    )
    existing_mols = []
    data = pd.read_csv(assay_path, index_col=0)
    assay_type = data.iloc[0]["activity.standard_type"]
    for index, row in data.iterrows():
        existing_mols.append(row["chembl_ligand.chembl_id"])
    # print(assay_path.stem[:-2])
    maybe_new_molecules = get_activities_smiles(assay_path.stem[:-2], assay_type)
    # print(maybe_new_molecules)
    data = get_components_for_target("CHEMBL3038498")
    print(data)
    sequences = get_target_clusters(["CHEMBL3038498"])
    print(sequences)


# ====================================================================#

if __name__ == "__main__":
    main()

# ====================================================================#


"""
Potential missing data - darf leben
Potential author error - darf leben
Manually validated - ok
Potential transcription error - darf leben
Outside typical range - killen
Non standard unit for type - killen
Author confirmed error - killen
"""
