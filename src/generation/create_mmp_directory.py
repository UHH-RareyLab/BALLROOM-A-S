import json
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from src.generation.parse_mmp_output import MMPOutputHandler
from src.utils.definitions import rdkit_path


def make_mmp_benchmark(
    start_path,
    end_path,
    rdkit=True,
    mmpdb=True,
    rdkit_hac=10,
    rdkit_ratio=0.3,
    min_core_hac=8,
):
    ### Here we copy the required data from the existing benchmark

    if end_path.exists():
        shutil.rmtree(end_path.as_posix())

    if not end_path.exists():
        end_path.mkdir()

    for target in start_path.glob("*"):
        if target.is_file():
            continue
        intermediate_path = end_path / target.name
        if not intermediate_path.exists():
            intermediate_path.mkdir()
        for assay in target.glob("*.csv"):
            original_smiles_path = assay.parent / (assay.stem + ".smi")
            original_data_path = assay.parent / (assay.stem + ".csv")
            extra_smiles_path = assay.parent / (assay.stem + "_additional.smiles")
            extra_data_path = assay.parent / (assay.stem + "_additional.tsv")
            original_smiles = pd.read_csv(original_smiles_path, sep=" ", header=None)
            original_data = pd.read_csv(
                original_data_path, dtype=str, sep=",", index_col=0
            )
            original_data["naomi_name"] = (
                original_data["pdb"] + "_" + original_data["name"]
            )
            activity_data = original_data[
                [
                    "pdb_ligand.naomi_unique_smiles_chiral",
                    "activity.standard_value",
                    "naomi_name",
                    "activity.standard_type",
                    "activity.standard_units",
                    "activity.chembl_activity_id",
                ]
            ]
            activity_data.columns = [
                "SMILES",
                "Activity",
                "Name",
                "Activity Type",
                "Activity Unit",
                "Activity_id",
            ]
            type = activity_data["Activity Type"].iloc[0]
            unit = activity_data["Activity Unit"].iloc[0]
            if extra_smiles_path.exists():
                extra_smiles = pd.read_csv(extra_smiles_path, sep=" ", header=None)
                extra_data = pd.read_csv(
                    extra_data_path, sep="\t", dtype=str, index_col=0
                )
                extra_data.columns = ["SMILES", "Activity", "Name", "Activity_id"]
                extra_data["Activity Type"] = type
                extra_data["Activity Unit"] = unit
                original_smiles = pd.concat([original_smiles, extra_smiles])
                activity_data = pd.concat([activity_data, extra_data])
            new_path = intermediate_path / (assay.stem + ".smi")
            new_path_data = intermediate_path / (assay.stem + ".csv")
            # Sometimes, there are multiple values for the same molecule in the data even though same endpoints are given
            # In that case, I want to discard the whole assay because I dont trust it
            # I figure that out by deduplicating the smiles and if one is double discard the whole assay
            deduplicated_smiles = original_smiles.drop_duplicates(ignore_index=True)
            if deduplicated_smiles.shape[0] != original_smiles.shape[0]:
                print(
                    f"Assay {assay.name} has a problem with multiple activities to the same molecule, discarding...."
                )
                continue
            original_smiles.to_csv(
                new_path.as_posix(), header=False, index=False, sep=" "
            )
            activity_data.to_csv(new_path_data.as_posix(), index=False)

    ### Generate the actual MMP data on the molecules
    base_rdkit_path = rdkit_path / "Contrib" / "mmpa"
    fragment_path = base_rdkit_path / "rfrag.py"
    indexing_path = base_rdkit_path / "indexing.py"
    indexing_path = Path(
        "/local/gutermuth/backup/bob_analyse_skripte/gutermuth/own_rdkit/indexing.py"
    )
    mmpdb_path = Path("/work/gutermuth/mmpdb-2.1/mmpdb")
    print(base_rdkit_path)
    print(fragment_path.exists(), indexing_path.exists())
    python_path = os.path.abspath("/scratch/gutermuth/anaconda3/envs/bigboy/bin/python")

    for target in end_path.glob("*"):
        if target.is_file():
            continue
        for assay in target.glob("*.smi"):
            if not mmpdb:
                print("rdkit non mmpdb")
                fragment_command = [
                    "python",
                    str(fragment_path),
                    "<",
                    str(assay),
                    ">",
                    str(assay.parent / (assay.stem + ".frag")),
                ]
                fragment_command = " ".join(fragment_command)
                subprocess.run(fragment_command, shell=True)
                index_command_base = [
                    "python",
                    str(indexing_path),
                    "<",
                    str(assay.parent / (assay.stem + ".frag")),
                    ">",
                    str(assay.parent / (assay.stem + "_heavy_atoms.mmp")),
                ]
                # subprocess.run(" ".join(index_command_base), shell=True)
                index_command_ratio = [
                    "python",
                    str(indexing_path),
                    "<",
                    str(assay.parent / (assay.stem + ".frag")),
                    "-r",
                    "0.3",
                    ">",
                    str(assay.parent / (assay.stem + "_ratio.mmp")),
                ]
                subprocess.run(" ".join(index_command_ratio), shell=True)
            else:
                print(assay)
                print("rdkit mmpdb")
                fragment_command = [
                    python_path,
                    mmpdb_path,
                    "fragment",
                    assay,
                    "-o",
                    str(assay.parent / (assay.stem + ".frag")),
                ]
                subprocess.run(fragment_command, check=True)
                index_command = [
                    python_path,
                    mmpdb_path,
                    "index",
                    str(assay.parent / (assay.stem + ".frag")),
                    "--max-variable-ratio",
                    str(rdkit_ratio),
                    "--max-variable-heavies",
                    str(rdkit_hac),
                    "-o",
                    str(assay.parent / (assay.stem + "_ratio.mmp")),
                    "--out",
                    "csv",
                ]
                subprocess.run(index_command, check=True)
    print("Finished calculations")
    for target in end_path.glob("*"):
        if target.is_file():
            continue
        for mmp_data in target.glob("*.mmp"):
            current = MMPOutputHandler(
                mmp_data, rdkit=rdkit, mmpdb=mmpdb, minimum_core=min_core_hac
            )
            current.run()

    print("Finished mmp finding")
    ### Generate CSVs with data in MMP directories
    for target in end_path.glob("*"):
        if target.is_file():
            continue
        for directory in target.glob("*"):
            if directory.is_file():
                continue
            for mmp_run in directory.glob("*.smi"):
                mmp_smiles_data = pd.read_csv(
                    mmp_run, sep=" ", header=None, index_col=None
                )
                assay_data_path = target / (directory.name.split("_")[0] + ".csv")
                assay_data = pd.read_csv(assay_data_path)
                activities = []
                activity_types = []
                activity_units = []
                activity_ids = []
                for index, row in mmp_smiles_data.iterrows():
                    name = row[1]
                    found = False
                    for index2, row2 in assay_data.iterrows():
                        name_assay = row2["Name"]
                        activity = row2["Activity"]
                        activity_type = row2["Activity Type"]
                        activity_unit = row2["Activity Unit"]
                        activity_id = row2["Activity_id"]
                        if name == name_assay:
                            activities.append(activity)
                            activity_units.append(activity_unit)
                            activity_types.append(activity_type)
                            activity_ids.append(activity_id)
                            found = True
                            break
                    if not found:
                        raise RuntimeError
                mmp_smiles_data["Activity"] = activities
                mmp_smiles_data["Activity Type"] = activity_types
                mmp_smiles_data["Activity Unit"] = activity_units
                mmp_smiles_data["Activity_id"] = activity_ids
                mmp_smiles_data.columns = [
                    "SMILES",
                    "Name",
                    "Activity",
                    "Activity Type",
                    "Activity Unit",
                    "Activity_id",
                ]
                min_activity = 10000
                max_activity = -1000
                for index, row in mmp_smiles_data.iterrows():
                    activity = row["Activity"]
                    if activity < min_activity:
                        min_activity = activity
                    if activity > max_activity:
                        max_activity = activity
                # if (orderOfMagnitude(max_activity) - orderOfMagnitude(min_activity)) >= 3:
                mmp_smiles_data.to_csv(
                    mmp_run.parent / (mmp_run.stem + ".csv"), index=None
                )
                # else:
                #    mmp_run.unlink()

    print("Finished generating mmp data ")
    ### Cleanup and results
    for target in end_path.glob("*"):
        for frag in target.glob("*.frag"):
            frag.unlink()
        for mmp in target.glob("*.mmp"):
            mmp.unlink()
        if target.is_file():
            continue
        assay_names = []
        number_dirs = len([x for x in target.glob("*") if x.is_dir()])
        if number_dirs == 0:
            shutil.rmtree(target)
            continue
        for directory in target.glob("*"):
            if directory.is_file():
                continue
            assay_names.append(directory.name.split("_")[0])
        for file in target.glob("*"):
            if file.is_dir():
                continue
            if not file.stem in assay_names:
                file.unlink()

    n_targets = 0
    n_assays = 0
    n_mmps = 0
    n_datapoints = 0
    for target in end_path.glob("*"):
        n_targets += 1
        current_assays = len([x for x in target.glob("*.csv")])
        n_assays += current_assays
        for assay in target.glob("*"):
            if assay.is_file():
                continue
            for mmp_data in assay.glob("*.csv"):
                n_mmps += 1
                full_data = pd.read_csv(mmp_data)
                n_datapoints += full_data.shape[0]

    print(n_targets, n_assays, n_mmps, n_datapoints)

    data_dict = {
        "n_targets": n_targets,
        "n_assays": n_assays,
        "n_mms": n_mmps,
        "n_datapoints": n_datapoints,
    }
    settings_path = end_path / "results.json"
    with open(settings_path, "w") as f:
        json.dump(data_dict, f)


def main():
    start_path = Path(
        "/local/gutermuth/backup/bob_analyse_skripte/gutermuth/benchmark_directory/mtz_certain_noncontradicting_elaborate/pair_benchmark_starting_point_strict_mols_merged"
    )
    end_path = Path(
        "/local/gutermuth/backup/bob_analyse_skripte/gutermuth/benchmark_directory/testaround_hussein_rea"
    )
    make_mmp_benchmark(
        start_path=start_path, end_path=end_path, rdkit=True, mmpdb=False
    )


if __name__ == "__main__":
    main()
