import json
from pathlib import Path

import pandas as pd

import src.utils.definitions as definitions
from src.utils.helper_functions import update_dict


def extract_jamda_dockings():
    dockings_csv_path = definitions.docking_output_path / "docking_df_jamda.csv"
    docking_data = pd.read_csv(dockings_csv_path, index_col=None, header=None)
    jamda_score_path = definitions.score_path / "jamda_scores.json"
    jamda_rmsd_path = definitions.score_path / "jamda_rmsds.json"
    docking_data_dict = {}
    rmsd_data_dict = {}
    for index, row in docking_data.iterrows():
        current_dict = {}
        current_docking = Path(row[3])
        pdb = current_docking.parent.stem.split("_")[0]
        ligand = current_docking.parent.stem[5:]
        csvs = [
            x
            for x in current_docking.glob("*.csv")
            if not x.stem.endswith("protonation")
        ]
        for csv in csvs:
            data = pd.read_csv(csv, sep=";")
            if data.shape[0] > 0:
                name = data.loc[0]["molecule_name"]
                score = data.loc[0]["score"]
                rmsd = data.loc[0]["rmsd"]
                if not pd.isna(rmsd):
                    rmsd_data_dict[name] = rmsd
                current_dict[name] = score
        docking_data_dict[pdb + "-" + ligand] = current_dict
    if jamda_score_path.exists():
        with open(jamda_score_path, "r") as f:
            current_jamda_data = json.load(f)
        complete_jamda_data, changes = update_dict(
            current_jamda_data, docking_data_dict
        )
        print(f"There have been {changes} updates to the jamda dict")
    else:
        complete_jamda_data = docking_data_dict
    with open(jamda_score_path, "w") as f:
        json.dump(complete_jamda_data, f)
    if jamda_rmsd_path.exists():
        with open(jamda_rmsd_path, "r") as f:
            current_jamda_data = json.load(f)
        complete_jamda_data, changes = update_dict(current_jamda_data, rmsd_data_dict)
        print(f"There have been {changes} updates to the jamda RMSD dict")
    else:
        complete_jamda_data = rmsd_data_dict
    with open(jamda_rmsd_path, "w") as f:
        json.dump(complete_jamda_data, f)


if __name__ == "__main__":
    extract_jamda_dockings()
