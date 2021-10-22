import shutil
import subprocess

import pandas as pd

import src.utils.definitions as definitions


def copy_and_extract_pdb(output_directory, pdb, ligandname):
    input_pdb_name = definitions.mirror_path / ("pdb" + pdb.lower() + ".ent.gz")
    output_pdb_name = output_directory / (pdb + ".pdb.gz")
    final_pdb_name = output_directory / (pdb + ".pdb")
    if not final_pdb_name.exists():
        shutil.copyfile(input_pdb_name, output_pdb_name)
        subprocess.run(["gzip", "-d", output_pdb_name])
    output_ligand_name = output_directory / (pdb + "_" + ligandname + "_ligand.sdf")
    task = [
        str(definitions.unicon_path),
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


def get_all_receptors_from_string(docking_string):
    final_list = []
    for receptor in docking_string.split(";"):
        if receptor != "":
            final_list.append(receptor)
    return final_list


def figure_out_dockings():
    base_output_path = definitions.docking_path
    base_output_path.mkdir(exist_ok=True)
    definitive_docking_dict = {}

    for benchmark in definitions.benchmarks:
        data_path = benchmark / "complete_data.csv"
        data = pd.read_csv(data_path, index_col=0)

        for index, row in data.iterrows():
            smiles1 = row["Smiles1"]
            name1 = row["ID1"]
            smiles2 = row["Smiles2"]
            name2 = row["ID2"]
            dockinto1 = get_all_receptors_from_string(row["Dockinto1"])
            dockinto2 = get_all_receptors_from_string(row["Dockinto2"])
            for receptor in dockinto1:
                if receptor in definitive_docking_dict.keys():
                    definitive_docking_dict[receptor].add(smiles1 + " " + name1)
                else:
                    definitive_docking_dict[receptor] = set()
                    definitive_docking_dict[receptor].add(smiles1 + " " + name1)
            for receptor in dockinto2:
                if receptor in definitive_docking_dict.keys():
                    definitive_docking_dict[receptor].add(smiles2 + " " + name2)
                else:
                    definitive_docking_dict[receptor] = set()
                    definitive_docking_dict[receptor].add(smiles2 + " " + name2)
    print(definitive_docking_dict)
    definitions.docking_output_path.mkdir(exist_ok=True)
    docking_data = []
    # final docking file ist immer PDB_path active_site_path ligand_file_path output_folder_path
    ligands_docked = 0
    for receptor in definitive_docking_dict.keys():
        ligands_docked += len(definitive_docking_dict[receptor])
        print(
            receptor,
            len(definitive_docking_dict[receptor]),
            definitive_docking_dict[receptor],
        )
        pdb = receptor.split("_")[0]
        ligand = receptor[5:]
        print(pdb, ligand)
        working_dir = definitions.docking_output_path / receptor
        working_dir.mkdir(exist_ok=True)
        pdb_path, active_site_path = copy_and_extract_pdb(
            output_directory=working_dir, pdb=pdb, ligandname=ligand
        )
        ligand_path = working_dir / "to_dock.smi"
        with open(ligand_path, "w") as f:
            for to_dock in definitive_docking_dict[receptor]:
                f.write(to_dock + "\n")
        docking_output = working_dir / "results"
        if not docking_output.exists():
            docking_output.mkdir()
        docking_data.append(
            [str(pdb_path)[:-3], active_site_path, ligand_path, docking_output]
        )

    docking_data = pd.DataFrame(docking_data)
    docking_data.to_csv(
        definitions.docking_output_path / "docking_data.csv", header=None, index=None
    )
