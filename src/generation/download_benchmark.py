import itertools
import shutil
import subprocess
import urllib.request

import pandas as pd
import progressbar

from src.utils.definitions import pdb_download_url, unicon_path, mirror_path
from src.utils.helper_functions import (
    extract_pdb_ligand_from_name,
    get_all_mmp_data_in_benchmark,
)


def download_pdb(pdb, path_to_download):
    url_to_request = pdb_download_url + pdb.upper() + ".pdb"
    urllib.request.urlretrieve(url_to_request, path_to_download)


def copy_and_extract_pdb(output_directory, pdb, ligandname):
    input_pdb_name = mirror_path / ("pdb" + pdb.lower() + ".ent.gz")
    output_pdb_name = output_directory / (pdb + ".pdb.gz")
    final_pdb_name = output_directory / (pdb + ".pdb")
    if not final_pdb_name.exists():
        shutil.copyfile(input_pdb_name, output_pdb_name)
        subprocess.run(["gzip", "-d", output_pdb_name])
    output_ligand_name = output_directory / (pdb + "_" + ligandname + "_ligand.sdf")
    task = [
        str(unicon_path),
        "-i",
        str(final_pdb_name),
        "-o",
        str(output_ligand_name),
        "--extract",
        ligandname,
    ]
    if not output_ligand_name.exists():
        subprocess.run(task, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return [output_pdb_name, output_ligand_name]


def download_pair_benchmark(
    benchmark_path, unicon_path=unicon_path, mirror_path=mirror_path
):
    iterator_target = 0
    widget = [
        "Downloading/Copying benchmark: ",
        progressbar.BouncingBar(),
        progressbar.AdaptiveETA(),
    ]
    assays = [y / "pair_data.csv" for y in benchmark_path.glob("*") if y.is_dir()]
    for assay in progressbar.progressbar(assays, max_value=len(assays), widgets=widget):
        assay_data = pd.read_csv(assay, index_col=0, dtype=str)
        for index, row in assay_data.iterrows():
            id1 = row["ID1"]
            id2 = row["ID2"]
            result1 = extract_pdb_ligand_from_name(id1)
            result2 = extract_pdb_ligand_from_name(id2)
            if result1:
                pdb = result1[0]
                name = result1[1][5:]
                copy_and_extract_pdb(assay.parent, pdb, name)
            if result2:
                pdb = result2[0]
                name = result2[1][5:]
                copy_and_extract_pdb(assay.parent, pdb, name)
    iterator_target += 1


def download_mmp_pair_benchmark(
    benchmark_path, unicon_path=unicon_path, mirror_path=mirror_path
):
    iterator_target = 0
    widget = [
        "Downloading/Copying benchmark: ",
        progressbar.BouncingBar(),
        progressbar.AdaptiveETA(),
    ]
    mmps = get_all_mmp_data_in_benchmark(benchmark_path)
    for mmp in progressbar.progressbar(mmps, max_value=len(mmps), widgets=widget):
        assay_data = pd.read_csv(mmp, index_col=0, dtype=str)
        for index, row in assay_data.iterrows():
            id1 = row["Name"]
            result1 = extract_pdb_ligand_from_name(id1)
            if result1:
                pdb = result1[0]
                name = result1[1][5:]
                copy_and_extract_pdb(mmp.parent.parent, pdb, name)
    iterator_target += 1


def download_benchmark(
    benchmark_path, unicon_path=unicon_path, mirror_path=mirror_path
):
    """
    Downloads (or rather copies...) the pdbs and extracts the sdfiles using unicon
    @param benchmark_path: Path to the benchmark that shall be used
    @param unicon_path: Path to unicon executable
    @param mirror_path: Path to the pdb mirror that shall be used
    @return:
    """
    iterator_target = 0
    widget = [
        "Downloading/Copying benchmark: ",
        progressbar.BouncingBar(),
        progressbar.AdaptiveETA(),
    ]
    assays = list(
        itertools.chain(
            *[
                [y for y in x.glob("*.csv")]
                for x in benchmark_path.glob("*")
                if x.is_dir()
            ]
        )
    )
    for assay in progressbar.progressbar(assays, max_value=len(assays), widgets=widget):
        assay_data = pd.read_csv(assay, index_col=0, dtype=str)
        for index, row in assay_data.iterrows():
            name = row["name"]
            pdb = str(row["pdb.pdb_id"])
            input_pdb_name = mirror_path / ("pdb" + pdb.lower() + ".ent.gz")
            output_pdb_name = assay.parent / (pdb + ".pdb.gz")
            final_pdb_name = assay.parent / (pdb + ".pdb")
            if not final_pdb_name.exists():
                shutil.copyfile(input_pdb_name, output_pdb_name)
                subprocess.run(["gzip", "-d", output_pdb_name])
            output_ligand_name = assay.parent / (pdb + "_" + name + "_ligand.sdf")
            task = [
                str(unicon_path),
                "-i",
                str(final_pdb_name),
                "-o",
                str(output_ligand_name),
                "--extract",
                name,
            ]
            if not output_ligand_name.exists():
                subprocess.run(
                    task, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
    iterator_target += 1
