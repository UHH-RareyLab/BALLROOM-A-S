import copy

import numpy as np
import pandas as pd

from src.utils.db_interaction import query_chembl


def alter_matrix(matrix, index1, index2, value=0):
    matrix[index1][index2] = value
    matrix[index2][index1] = value


def check_contradictions_between_pairs(pairs1, pairs2):
    contradicting = False
    n_contradictions = 0
    for pair1 in pairs1:
        for pair2 in pairs2:
            if pair1.id1 == pair2.id1 and pair1.id2 == pair2.id2:
                if (
                    pair1.activity1 > pair1.activity2
                    and pair2.activity1 > pair2.activity2
                ) or (
                    pair1.activity1 < pair1.activity2
                    and pair2.activity1 < pair2.activity2
                ):
                    continue
                else:
                    contradicting = True
                    n_contradictions += 1
            elif pair1.id1 == pair2.id2 and pair1.id2 == pair2.id1:
                if (
                    pair1.activity1 > pair1.activity2
                    and pair2.activity1 < pair2.activity2
                ) or (
                    pair1.activity1 < pair1.activity2
                    and pair2.activity1 > pair2.activity2
                ):
                    continue
                else:
                    contradicting = True
                    n_contradictions += 1
    return contradicting, n_contradictions


def check_if_matrix_works(matrix, indexes_to_delete):
    short_matrix = copy.deepcopy(matrix)
    short_matrix = np.delete(short_matrix, indexes_to_delete, 0)
    short_matrix = np.delete(short_matrix, indexes_to_delete, 1)
    for i in range(len(short_matrix)):
        for j in range(len(short_matrix)):
            if i == j:
                continue
            if short_matrix[i][j] != 1:
                return False
    return True


def delete_until_happy(number_dict, adjacency_matrix, label_list):
    adjaceny_matrix_work = copy.deepcopy(adjacency_matrix)
    perc_dict = {x: number_dict[x][0] / number_dict[x][1] for x in number_dict.keys()}
    sorted_dict = sorted(perc_dict.items(), key=lambda x: x[1], reverse=True)
    index = 0
    indexes_to_delete = []
    while True:
        if index >= len(sorted_dict):
            print("alles rausschmeissen?")
            return [x[0] for x in sorted_dict]
        if check_if_matrix_works(adjaceny_matrix_work, indexes_to_delete):
            return indexes_to_delete
        indexes_to_delete.append(label_list.index(sorted_dict[index][0]))
        index += 1


def check_targets_assays(benchmark_path, conflicting_assays):
    n_targets = 0
    n_assays = 0
    n_datapoints = 0
    n_pockets = 0
    n_ligands = 0
    n_pairs = 0
    data = pd.read_csv(benchmark_path / "complete_data.csv")
    all_targets = set()
    all_assays = set()
    all_data = set()
    all_ligands = set()
    all_pockets = set()
    for index, row in data.iterrows():
        target = row["Target"]
        assay = row["Assay"]
        if assay in conflicting_assays:
            continue
        n_pairs += 1
        id1 = row["ID1"]
        all_ligands.add(id1)
        id2 = row["ID2"]
        all_ligands.add(id2)
        dockinto1 = row["Dockinto1"]
        dockinto2 = row["Dockinto1"]
        datapoint1 = target + "_" + assay + "_" + id1
        datapoint2 = target + "_" + assay + "_" + id2
        if not target in all_targets:
            all_targets.add(target)
            n_targets += 1
        if not assay in all_assays:
            all_assays.add(assay)
            n_assays += 1
        if not datapoint1 in all_data:
            all_data.add(datapoint1)
            n_datapoints += 1
        if not datapoint2 in all_data:
            all_data.add(datapoint2)
            n_datapoints += 1
        if type(dockinto1) != str or type(dockinto2) != str:
            print("skipped")
            continue
        for receptor in dockinto1.split(";"):
            if not receptor == "":
                all_pockets.add(receptor)
        for receptor in dockinto2.split(";"):
            if not receptor == "":
                all_pockets.add(receptor)
    return (
        n_targets,
        n_assays,
        n_datapoints,
        len(all_pockets),
        len(all_ligands),
        n_pairs,
    )


def get_noncontradicting_set(pair_data, target_name, most_assays=True):
    label_assays = []
    for index, row in pair_data.iterrows():
        assay = row["Assay"]
        if assay not in label_assays:
            label_assays.append(assay)
    query = f"""SELECT chembl_id, variant_id FROM assays WHERE assays.chembl_id IN ({','.join([f"'{id}'" for id in label_assays])})"""
    variant_data = query_chembl(query)
    variant_assays = []
    for index, row in variant_data.iterrows():
        assay = row["chembl_id"]
        variant_id = row["variant_id"]
        if not pd.isna(variant_id):
            variant_assays.append(assay)
    label_assays = [x for x in label_assays if x not in variant_assays]
    size_matrix = len(label_assays)
    start_matrix = np.full((size_matrix, size_matrix), 1)
    n_pairs = 0
    n_contradictions = 0
    pair_data_target = {}
    contradicting_assays = set()
    number_dict = {}
    for index, row in pair_data.iterrows():
        n_pairs += 1
        id1 = row["ID1"]
        id2 = row["ID2"]
        activity1 = row["Activity1"]
        activity2 = row["Activity2"]
        assay = row["Assay"]
        if assay in variant_assays:
            continue
        if not assay in number_dict:
            number_dict[assay] = [0, 1]
        else:
            number_dict[assay][1] += 1
        if activity1 >= activity2:
            safe_string = id1
        else:
            safe_string = id2
        id_list = [id1, id2]
        id_list.sort()
        id_set = frozenset(id_list)
        if id_set in pair_data_target.keys():
            if pair_data_target[id_set][0] != safe_string:
                print(
                    f"Contradiction : {id1},{id2},{assay},{target_name},{safe_string},{pair_data_target[id_set]}"
                )
                n_contradictions += 1
                contradicting_assays.add(assay)
                contradicting_assays.add(pair_data_target[id_set][3])
                number_dict[assay][0] += 1
                number_dict[pair_data_target[id_set][3]][0] += 1
                alter_matrix(
                    start_matrix,
                    label_assays.index(assay),
                    label_assays.index(pair_data_target[id_set][3]),
                )
                pair_data_target[id_set][0] = "both wrong"
        else:
            pair_data_target[id_set] = [safe_string, id1, id2, assay, target_name]
    to_delete = delete_until_happy(number_dict, start_matrix, label_assays)
    assays_to_delete = [
        label_assays[i] for i in range(len(label_assays)) if i in to_delete
    ]
    return (
        assays_to_delete,
        n_pairs,
        n_contradictions,
        len(contradicting_assays),
        label_assays,
        start_matrix,
    )
