import json

from rdkit import Chem

import src.utils.definitions as definitions
from src.utils.helper_functions import update_dict


def extract_jamda_dockings():
    jamda_score_path = definitions.score_path / "jamda_scores_rescoring.json"
    jamda_rmsd_path = definitions.score_path / "jamda_rmsds_rescoring.json"
    rmsd_dict = {}
    score_dict = {}
    for complex in definitions.docking_output_path.glob("*"):
        if complex.is_file():
            continue
        pdb = complex.stem.split("_")[0]
        sdf_path = complex / "results" / "opt.sdf"
        mols = [x for x in Chem.SDMolSupplier(str(sdf_path))]
        assert len(mols) == 1
        score = mols[0].GetDoubleProp("SCORE")
        rmsd = mols[0].GetDoubleProp("RMSD")
        name = mols[0].GetProp("_Name")
        new_identifier = pdb + "_" + name
        first_identifier = pdb + "-" + name
        assert new_identifier == complex.name
        print(score, rmsd, name, complex.name, pdb)
        rmsd_dict[new_identifier] = rmsd
        score_dict[first_identifier] = {new_identifier: score}
    if jamda_score_path.exists():
        with open(jamda_score_path, "r") as f:
            current_jamda_data = json.load(f)
        complete_jamda_data, changes = update_dict(current_jamda_data, score_dict)
        print(f"There have been {changes} updates to the jamda dict")
    else:
        complete_jamda_data = score_dict
    with open(jamda_score_path, "w") as f:
        json.dump(complete_jamda_data, f)
    if jamda_rmsd_path.exists():
        with open(jamda_rmsd_path, "r") as f:
            current_jamda_data = json.load(f)
        complete_jamda_data, changes = update_dict(current_jamda_data, rmsd_dict)
        print(f"There have been {changes} updates to the jamda RMSD dict")
    else:
        complete_jamda_data = rmsd_dict
    with open(jamda_rmsd_path, "w") as f:
        json.dump(complete_jamda_data, f)
    return


if __name__ == "__main__":
    extract_jamda_dockings()
