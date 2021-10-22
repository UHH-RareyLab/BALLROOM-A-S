import json

import src.utils.definitions as definitions
from src.utils.helper_functions import update_dict


def get_names_from_to_dock(to_dock_path):
    with open(to_dock_path, "r") as f:
        lines = f.readlines()
    names = []
    iterator = 0
    for line in lines:
        iterator += 1
        line = line.strip()
        names.append([line.split(" ")[1], iterator])
    return names


def get_score_from_output_file(output_file_path):
    with open(output_file_path, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith("REMARK VINA RESULT:"):
                result = line.split()[3]
                return result
    return 10000


def get_ligand_dict_from_file(to_dock_file):
    counter = 1
    ligand_dict = {}
    with open(to_dock_file, "r") as f:
        for line in f.readlines():
            line = line.strip()
            ligand_name = line.split()[1]
            ligand_dict[counter] = ligand_name
            counter += 1
    return ligand_dict


mmp_scores_vina = {}

benchmark_path = definitions.docking_output_path
# benchmark_path = Path("/scratch/gutermuth/all_benchmark_dockings/presubmission_activityfinder_bugfix") / "general_dockings"
print(benchmark_path)

vina_scores = {}
for target in benchmark_path.glob("*"):
    if target.is_file():
        continue
    to_dock_file = target / "to_dock.smi"
    ligand_dict = get_ligand_dict_from_file(to_dock_file)
    for ligandfile in target.glob("*.sdf"):
        ligand_name_split = ligandfile.name.split("_")
        pdb = ligand_name_split[0]
        complex = ligandfile.name[:-11]
        pdb = complex.split("_")[0]
        ligand = complex[5:]
        complex = pdb + "-" + ligand
        output_folder = target / f"{pdb}_vina_prep" / "vina_output"
        for file in output_folder.glob("*"):
            score = get_score_from_output_file(file)
            counter = int(file.stem.split("_")[-1])
            print(score, counter, ligand_dict[counter])
            if not complex in vina_scores:
                vina_scores[complex] = {}
            vina_scores[complex][ligand_dict[counter]] = score


print(sum(1 for x in vina_scores.keys() if vina_scores[x] == 10000))
print(sum(1 for x in vina_scores.keys() if vina_scores[x] != 10000))
vina_score_path = definitions.score_path / "Vina.json"

if vina_score_path.exists():
    with open(vina_score_path, "r") as f:
        current_vina_data = json.load(f)
    complete_vina_data, changes = update_dict(current_vina_data, vina_scores)
    print(f"There have been {changes} updates to the Vina dict")
else:
    complete_vina_data = vina_scores
with open(vina_score_path, "w") as f:
    json.dump(complete_vina_data, f)
