import itertools
import shutil
import subprocess

import numpy as np
import pandas as pd
import progressbar
from scipy.sparse.csgraph import connected_components

from src.utils.helper_functions import extract_pdb_ligand_from_name


def get_stuff_from_label(label):
    """
    Extracts the pdb and ligand from a label of format [pdb]_[ligand]_ligand
    @param label: label to extract
    @return: list of pdb and ligand
    """
    pdb = label[0:4]
    ligand = label[5:-7]
    return [pdb, ligand]


def get_index_of_label(pdb, ligand, label_list):
    """
    Retrieves the index in a label_list of a pdb_ligand pair
    @param pdb: pdb code of interest
    @param ligand: ligand code of interest
    @param label_list: label list
    @return: index in label list or None
    """
    search_name = pdb + "_" + ligand + "_ligand"
    for index, value in enumerate(label_list):
        if search_name == value:
            return index
    return None


def initiate_adjacency_matrix(n):
    """
    Initiates an symetric adjacency matrix of size size
    @param n: size of adjacency matrix
    @return: initiated adjacency matrix
    """
    matrix = np.full((n, n), 0)
    for i in range(n):
        matrix[i][i] = 1
    return matrix


def create_adjacency_matrix(assay_subdirectory):
    """
    Creates and updates a adjacency matrix of pocket indexes in a given assay subdirectory
    @param assay_subdirectory: assoy of interest
    @return: list of the adjacency matrix and the labels of said matrix
    """
    runs = [x for x in (assay_subdirectory / "siena_runs").glob("*")]
    labels = [x.name for x in runs]

    a_matrix = initiate_adjacency_matrix(len(runs))
    for run in runs:
        search_pdb, search_ligand = get_stuff_from_label(run.name)
        index_search = get_index_of_label(search_pdb, search_ligand, labels)
        search_index_update_matrix(
            run / "resultStatistic.csv", labels, a_matrix, index_search
        )
    return [a_matrix, labels]


def calc_connected_components_assay(assay_subdirectory):
    """
    Calculates the connected components of the pocket information of an assay
    @param assay_subdirectory: assay of interest
    @return: list with the number of components, their order and the respective labels
    """
    matrix, labels = create_adjacency_matrix(assay_subdirectory=assay_subdirectory)
    number_components, order_components = connected_components(matrix)
    return [number_components, order_components, labels]


def check_and_make_dir(path):
    """
    Helper function to check if a directory exists and if not, create it
    @param path: directory to check
    @return: Nothing
    """
    if not path.exists():
        path.mkdir()
    else:
        shutil.rmtree(path)
        path.mkdir()


def extract_pdbs_ligands(assay_path):
    """
    Extracts all pdbs and ligands from an assay csv
    @param assay_path: Path to the csv of an assay
    @return: list with sets of all pdbs and ligands
    """
    pdbs = set()
    ligands = set()
    data = pd.read_csv(assay_path, index_col=0, dtype=str)
    for index, row in data.iterrows():
        pdb = row["pdb"] + ".pdb"
        ligand = row["pdb"] + "_" + row["name"] + "_ligand" + ".sdf"
        pdbs.add(pdb)
        ligands.add(ligand)
    return [pdbs, ligands]


def extract_pdbs_ligands_pairs(pair_data_path):
    pdbs = set()
    ligands = set()
    print(pair_data_path)
    data = pd.read_csv(pair_data_path, dtype=str, sep=",")
    for index, row in data.iterrows():
        print(row)
        dockinto1 = row["Dockinto1"]
        dockinto2 = row["Dockinto2"]
        print(dockinto1, dockinto2)
        for id1 in dockinto1.split(";"):
            result1 = extract_pdb_ligand_from_name(id1)
            if result1:
                pdbs.add(result1[0] + ".pdb")
                ligands.add(result1[1] + "_ligand.sdf")
        for id2 in dockinto2.split(";"):
            result2 = extract_pdb_ligand_from_name(id2)
            if result2:
                pdbs.add(result2[0] + ".pdb")
                ligands.add(result2[1] + "_ligand.sdf")
    return [pdbs, ligands]


def extract_pdbs_ligands_pairs_mmps(pair_data_path):
    pdbs = set()
    ligands = set()
    data = pd.read_csv(pair_data_path, dtype=str)
    for index, row in data.iterrows():
        id1 = row["Name"]
        result1 = extract_pdb_ligand_from_name(id1)
        if result1:
            pdbs.add(result1[0] + ".pdb")
            ligands.add(result1[1] + "_ligand.sdf")
    return [pdbs, ligands]


def make_siena_clusters(
    path_benchmark, path_siena, path_sienatools, benchmark="normal"
):
    """
    Function that creates the Siena clusters
    @param path_benchmark:
    @param path_siena:
    @param path_sienatools:
    @return:
    """
    complete_data = []
    widget = [
        "Creating Siena Clusters: ",
        progressbar.BouncingBar(),
        progressbar.AdaptiveETA(),
    ]
    if benchmark == "normal":
        assays = list(
            itertools.chain(
                *[
                    [y for y in x.glob("*.csv")]
                    for x in path_benchmark.glob("*")
                    if x.is_dir()
                ]
            )
        )
    elif benchmark == "pairs":
        assays = [y / "pair_data.csv" for y in path_benchmark.glob("*") if y.is_dir()]
    print(assays)
    # Iterate over all assays
    for assay in progressbar.progressbar(assays, max_value=len(assays), widgets=widget):
        if benchmark == "normal":
            pdbs, ligands = extract_pdbs_ligands(assay)
        elif benchmark == "pairs":
            pdbs, ligands = extract_pdbs_ligands_pairs(assay)
        chembl_assay_directory = assay.parent / assay.stem
        # Create folder structure
        check_and_make_dir(chembl_assay_directory)
        check_and_make_dir(chembl_assay_directory / "pdbs")
        check_and_make_dir(chembl_assay_directory / "siena_runs")
        for pdb in pdbs:
            startpdb = chembl_assay_directory.parent / pdb
            symlink = chembl_assay_directory / "pdbs" / pdb
            if not symlink.exists():
                (symlink).symlink_to(startpdb)
        # Create Database
        create_database = [
            str(path_sienatools),
            "-d",
            str(chembl_assay_directory / "pdbs"),
            "-b",
            str(chembl_assay_directory / "siena.db"),
            "-f",
            "1",
        ]
        subprocess.run(
            create_database, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        print(ligands)
        labels = [x[:-4] for x in ligands]
        print(assay.parent, labels, ligands, pdbs)
        # Query database and save results
        a_matrix = initiate_adjacency_matrix(len(ligands))
        for index, ligand in enumerate(ligands):
            search_pdb, search_ligand = get_stuff_from_label(ligand[:-4])
            print(search_pdb, search_ligand)
            index_search = get_index_of_label(search_pdb, search_ligand, labels)
            ligand_file = assay.parent / ligand
            pdb_file = assay.parent / (ligand[0:4] + ".pdb")
            output_directory = chembl_assay_directory / "siena_runs" / ligand[:-4]
            check_and_make_dir(output_directory)
            query_command = [
                str(path_siena),
                "-p",
                str(pdb_file),
                "-l",
                str(ligand_file),
                "-b",
                str(chembl_assay_directory / "siena.db"),
                "-o",
                str(output_directory),
                "-i",
                "0.7",
            ]
            if not (output_directory / "resultStatistic.csv").exists():
                subprocess.run(
                    query_command,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
            search_index_update_matrix(
                output_directory / "resultStatistic.csv", labels, a_matrix, index_search
            )
            # Remove unneccesary files to avoid file system bloat (sadly no option to do this in the CLI)
            if (output_directory / "ensemble").exists():
                shutil.rmtree(output_directory / "ensemble")
            if (output_directory / "ligands").exists():
                shutil.rmtree(output_directory / "ligands")
            # Calc important data
            n_pockets, indexes = connected_components(a_matrix)
            if n_pockets == 1:
                # stop the queries if all pockets are in the same cluster already
                break
        complete_data.append(
            [assay.parent.name, assay.stem, n_pockets, indexes, labels]
        )
    complete_data = pd.DataFrame(
        complete_data,
        columns=["Target", "Assay", "N_Pockets", "Pocket_Indexes", "Pocket_labels"],
    )
    complete_data.to_csv(path_benchmark / "pocket_data.csv")


def search_index_update_matrix(data_path, labels, a_matrix, index_search):
    """
    Searches the indexes of siena results and updates adjacency matrix
    @param data_path: Path to the results of a siena run
    @param labels: labels present in the siena run
    @param a_matrix: adjacency matrix
    @param index_search: The index with which we search (the query pocket)
    @return: Updates the adjacency matrix directly...
    """
    data = pd.read_csv(data_path, sep=";", dtype=str)
    for index, row in data.iterrows():
        if pd.isnull(row["Ligand PDB code"]):
            continue
        pdb = row["PDB code"]
        ligand = row["Ligand PDB code"]
        index_found = get_index_of_label(pdb, ligand, labels)
        if index_found == None:
            continue
        if index_search != index_found and a_matrix[index_search][index_found] == 0:
            a_matrix[index_search][index_found] = 1
            a_matrix[index_found][index_search] = 1
