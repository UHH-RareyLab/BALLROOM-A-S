import itertools

import pandas as pd
import progressbar


def finalize_siena_benchmark(benchmark_path, min_datapoints):
    """
    Uses the generated pocket data and filters the preliminary datasets present in the benchmark path according to it
    @param benchmark_path: path to the benchmark
    @return: filtered dataframes
    """
    data = pd.read_csv(benchmark_path / "pocket_data.csv", index_col=0, dtype=str)
    keepers = calc_keepers(data, min_datapoints)
    dataframes = {}
    assays = list(
        itertools.chain(
            *[
                [y for y in x.glob("*.csv")]
                for x in benchmark_path.glob("*")
                if x.is_dir()
            ]
        )
    )
    for assay in assays:
        key = assay.parent.name + "_" + assay.stem
        fitting_keys = [x for x in keepers.keys() if x.startswith(key)]
        for working_comb in fitting_keys:
            current_keepers = keepers[working_comb]
            to_keep = []
            raw_data = pd.read_csv(assay, index_col=0, dtype=str)
            for index, row in raw_data.iterrows():
                pdb = row["pdb"]
                ligand = row["name"]
                new_key = pdb + "_" + ligand + "_ligand"
                if new_key in current_keepers:
                    to_keep.append(index)
            filtered_dataframe = raw_data.loc[to_keep]
            key = (assay.parent.name, assay.stem, working_comb.split(".")[1])
            dataframes[key] = filtered_dataframe
    return dataframes


def most_frequent(list):
    """
    Finds the most frequent entry in a list
    @param list: list of interest
    @return: most frequent entry
    """
    return max(set(list), key=list.count)


def of_interest(test_list, min_obs):
    """
    Extracts all elements from a list that appear a minimum number of times
    @param test_list: List to extract elements from
    @param min_obs: minimal number of observations
    @return: dictionary of the kept elements in the list
    """
    d = {}
    for x in test_list:
        if x in d.keys():
            d[x] += 1
        else:
            d[x] = 1
    kept = {x: d[x] for x in d.keys() if d[x] >= min_obs}
    return kept


def calc_keepers(data, min_observations):
    """
    Calculates the pockets to keep in a siena run
    @param data: Siena data to check
    @param min_observations: minimum number of observations to find
    @return: dictionary with data to keep
    """
    n_kept = 0
    n_kept_true = 0
    target_assay_keepers = {}
    n_thrownaway = 0
    for index, row in data.iterrows():
        target = row["Target"]
        assay = row["Assay"]
        pocket_indexes = [
            int(x)
            for x in row["Pocket_Indexes"].replace("[", "").replace("]", "").split()
        ]
        pocket_labels = [
            x
            for x in row["Pocket_labels"]
            .replace("[", "")
            .replace("]", "")
            .replace("'", "")
            .replace(",", "")
            .split()
        ]
        all_interest = of_interest(pocket_indexes, min_observations)
        keepers = {x: [] for x in all_interest.keys()}
        for i, current_index in enumerate(pocket_indexes):
            if current_index in all_interest.keys():
                keepers[current_index].append(pocket_labels[i])
                n_kept += 1
            else:
                n_thrownaway += 1
        for final_index, key in enumerate(keepers.keys()):
            if len(keepers[key]) > min_observations:
                final_key = target + "_" + assay + "." + str(final_index)
                target_assay_keepers[final_key] = keepers[key]
                n_kept_true += len(keepers)
    return target_assay_keepers


def calc_active_site_keys(dataframe_dict):
    """
    Calculates all relevant active site keys for later extraction of data
    @param dataframe_dict: dictionary of dataframes
    @return: key map of active site keys initialized to 0 (no mutations)
    """
    key_map = {}
    for dataset in dataframe_dict:
        for index, row in dataframe_dict[dataset].iterrows():
            pdb = row["pdb.pdb_id"]
            ligand = row["name"]
            target = row["target.chembl_id"]
            assay = row["assay.chembl_id"]
            blast_id = row["blast_match.bmid"]
            if blast_id[-2] == ".":
                blast_id = blast_id[:-2]
            final_key = pdb + "-" + ligand + "-" + target + "-" + assay + "-" + blast_id
            if final_key in key_map.keys():
                continue
            else:
                key_map[final_key] = 0
    return key_map


def refine_active_site_data(raw_active_site_data_path, output_path):
    refined_active_site_data = {}
    length = sum(1 for _ in open(raw_active_site_data_path, "r"))
    chunksize = 10000
    steps = int(length / chunksize)
    widget = [
        "Extracting mutation data: ",
        progressbar.BouncingBar(),
        progressbar.AdaptiveETA(),
    ]
    for chunk in progressbar.progressbar(
        pd.read_csv(
            raw_active_site_data_path, index_col=0, dtype=str, chunksize=chunksize
        ),
        max_value=steps,
        widgets=widget,
    ):
        for index, row in chunk.iterrows():
            pdb = row["pdb.pdb_id"]
            ligand = row["active_site.naomi_ligand_id"]
            target = row["target.chembl_id"]
            assay = row["assay.chembl_id"]
            blast_id = row["blast_match.bmid"]
            final_key = pdb + "-" + ligand + "-" + target + "-" + assay + "-" + blast_id
            if final_key in refined_active_site_data.keys():
                refined_active_site_data[final_key] += 1
            else:
                refined_active_site_data[final_key] = 1
    with open(output_path, "w") as f:
        f.write("PDB, Ligand, Target, Assay, blast_id, N_mutations\n")
        for key in refined_active_site_data.keys():
            current_list = [x for x in key.split("-")]
            current_list.append(refined_active_site_data[key])
            print(current_list)
            f.write(
                f"{current_list[0]},{current_list[1]},{current_list[2]},{current_list[3]},{current_list[4]},{current_list[5]}\n"
            )
    return refined_active_site_data


def extract_mutation_data(key_map, mutation_data):
    """
    Extracts active_site_mutation_data
    @param key_map: the keys of active sites that are relevant
    @param active_site_data_path: path to the active site data
    @return: updates the key map with mutation numbers
    """
    for real_key in key_map.keys():
        if real_key in mutation_data.keys():
            key_map[real_key] = mutation_data[real_key]


def filter_dataframes_active_site(dataframes, max_mutations=0):
    """
    Filters dataframes according to if they have mutations in their active site
    @param dataframes: Dataframes to be filtered
    @param active_site_map: Map of active site mutation data
    @param max_mutations: maximum number of mutations allowed
    @return: Changes the dataframes inplace
    """
    for dataframe in dataframes:
        to_keep = []
        for index, row in dataframes[dataframe].iterrows():
            pdb = row["pdb.pdb_id"]
            ligand = row["name"]
            target = row["target.chembl_id"]
            assay = row["assay.chembl_id"]
            blast_id = row["blast_match.bmid"]
            if blast_id[-2] == ".":
                blast_id = blast_id[:-2]
            n_mut = row["N_mutations"]
            if int(n_mut) <= max_mutations:
                to_keep.append(index)
            else:
                print(
                    "We filtered something",
                    pdb,
                    ligand,
                    target,
                    assay,
                    blast_id,
                    n_mut,
                    dataframe,
                )
        dataframes[dataframe] = dataframes[dataframe].loc[to_keep]
