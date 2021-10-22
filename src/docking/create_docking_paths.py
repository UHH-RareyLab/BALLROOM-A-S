import json

import pandas as pd

import src.utils.definitions as definitions


def figure_out_vina_dockings():
    list_to_dock = []

    for target in definitions.docking_output_path.glob("*"):
        if target.is_file() or target.name == "logs":
            continue
        print(target, target.name)
        for prep_folder in target.glob("*_vina_prep"):
            if prep_folder.is_dir():
                # print(prep_folder)
                box_path = prep_folder / "box_dimensions.json"
                output_folder = prep_folder / "vina_output"
                pdb = [x for x in prep_folder.glob("*.pdbqt")][0]
                for ligand in (prep_folder / "prepped_ligands").glob("*.pdbqt"):
                    stuff_to_dock = [
                        pdb.as_posix(),
                        ligand.as_posix(),
                        box_path.as_posix(),
                        output_folder.as_posix(),
                    ]
                    print(stuff_to_dock)
                    list_to_dock.append(stuff_to_dock)
    docking_df = pd.DataFrame(list_to_dock)
    docking_df.to_csv(
        definitions.docking_output_path / "docking_df_vina.csv",
        index=None,
        header=None,
        sep=" ",
    )


def figure_out_dock_dockings():
    dock_dockings = []
    for target in definitions.docking_output_path.glob("*"):
        if target.is_file() or target.name == "logs":
            continue
        dock_dockings.append([target])

    docking_df = pd.DataFrame(dock_dockings)
    docking_df.to_csv(
        definitions.docking_output_path / "docking_df_dock.csv",
        index=None,
        header=None,
        sep=" ",
    )


def initialize_existing_scores():
    dock_scores = []
    vina_scores = []
    jamda_scores = []
    scoring_path = definitions.score_path
    with open(scoring_path / "DOCK.json") as f:
        dock_scores = json.load(f)
    with open(scoring_path / "Vina.json") as f:
        vina_scores = json.load(f)
    with open(scoring_path / "jamda_scores.json") as f:
        jamda_scores = json.load(f)
    return jamda_scores, vina_scores, dock_scores


def target_present_in_scores(target, scores, ligand_dict):
    final_target = (
        target.name.split("_")[0] + "-" + "_".join(target.name.split("_")[1:])
    )
    for ligand_id in ligand_dict:
        if not final_target in scores:
            return False
        if not ligand_dict[ligand_id] in scores[final_target]:
            return False
    return True


def check_vina_score_is_present(target, ligand, vina_scores, ligand_dict):
    true_target = target.name.split("_")[0] + "-" + "_".join(target.name.split("_")[1:])
    ligand_id = int(ligand.stem.split("_")[-1])
    if not true_target in vina_scores:
        return False
    if not ligand_dict[ligand_id] in vina_scores[true_target]:
        return False
    return True


def get_ligands_from_smi_file(target):
    file_path = target / "to_dockprepared.smi"
    if not file_path.exists():
        raise RuntimeError("No prepared smiles file found")
    ligand_data = pd.read_csv(
        file_path, sep=" ", header=None, index_col=False, names=["Smiles", "Name"]
    )
    ligand_dict = {}
    for index, row in ligand_data.iterrows():
        ligand_dict[index + 1] = row["Name"]
    return ligand_dict


def write_docking_dfs():
    vina_dockings = []
    dock_dockings = []
    jamda_dockings = []
    jamda_scores, vina_scores, dock_scores = initialize_existing_scores()
    for target in definitions.docking_output_path.glob("*"):
        print(f"analyzing target {target}")
        if target.is_file() or target.name == "logs":
            continue
        ligand_dict = get_ligands_from_smi_file(target)
        if not target_present_in_scores(target, dock_scores, ligand_dict):
            print(f"DOCK dockings need to be repeated for target {target}")
            dock_dockings.append([target])
        for prep_folder in target.glob("*_vina_prep"):
            if prep_folder.is_dir():
                box_path = prep_folder / "box_dimensions.json"
                output_folder = prep_folder / "vina_output"
                pdb = [x for x in prep_folder.glob("*.pdbqt")][0]
                for ligand in (prep_folder / "prepped_ligands").glob("*.pdbqt"):
                    to_dock_vina = [
                        pdb.as_posix(),
                        ligand.as_posix(),
                        box_path.as_posix(),
                        output_folder.as_posix(),
                    ]
                    if not check_vina_score_is_present(
                        target, ligand, vina_scores, ligand_dict
                    ):
                        print(
                            f"VINA docking needs to be repeated for target {target} and ligand {ligand}"
                        )
                        vina_dockings.append(to_dock_vina)
        for pdb in target.glob("*.pdb"):
            ligands = [x for x in target.glob("*.sdf") if x.name.startswith(pdb.stem)]
            assert len(ligands) == 1
            ligand = ligands[0]
            ligands_to_dock = target / "to_dockprepared.smi"
            results = target / "results"
            assert ligands_to_dock.exists()
            assert results.exists()
            if (
                not target_present_in_scores(target, jamda_scores, ligand_dict)
            ) or True:
                print(f"JAMDA docking needs to be repeated for target {target}")
                jamda_dockings.append(
                    [
                        pdb.as_posix(),
                        ligand.as_posix(),
                        ligands_to_dock.as_posix(),
                        results.as_posix(),
                    ]
                )
    docking_df = pd.DataFrame(dock_dockings)
    print(f"Dock dockings to repeat: {len(dock_dockings)}")
    docking_df.to_csv(
        definitions.docking_output_path / "docking_df_dock.csv",
        index=None,
        header=None,
        sep=" ",
    )
    docking_df = pd.DataFrame(vina_dockings)
    print(f"Vina dockings to repeat: {len(vina_dockings)}")
    docking_df.to_csv(
        definitions.docking_output_path / "docking_df_vina.csv",
        index=None,
        header=None,
        sep=" ",
    )
    docking_df = pd.DataFrame(jamda_dockings)
    print(f"JAMDA dockings to repeat: {len(jamda_dockings)}")
    docking_df.to_csv(
        definitions.docking_output_path / "docking_df_jamda.csv",
        index=None,
        header=None,
        sep=",",
    )
