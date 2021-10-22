import json
from pathlib import Path

import src.utils.definitions as definitions
from src.utils.helper_functions import update_dict

current_path = Path(
    "/work/gutermuth/boltz_tryaround/output_dir/boltz_results_newest_boltz_dockings/predictions"
)
benchmark_data_dict_affinity = {}
benchmark_data_dict_probability = {}

for folder in current_path.glob("*"):
    if folder.is_file():
        continue
    result_path = folder
    affinity_results = [x for x in result_path.glob("affinity*.json")]
    if len(affinity_results) != 1:
        print(
            f"Folder {folder.name} skipped because there are {len(affinity_results)} results"
        )
        continue
    pdb = affinity_results[0].name.split("_")[1]
    molecule = affinity_results[0].name[14:-5]
    with open(str(affinity_results[0]), "r") as f:
        data = json.load(f)
    current_affinity = data["affinity_pred_value"]
    current_probability = data["affinity_probability_binary"]
    if pdb in benchmark_data_dict_affinity.keys():
        benchmark_data_dict_probability[pdb][molecule] = current_probability
        benchmark_data_dict_affinity[pdb][molecule] = current_affinity
    else:
        benchmark_data_dict_probability[pdb] = {molecule: current_probability}
        benchmark_data_dict_affinity[pdb] = {molecule: current_affinity}
print(len(benchmark_data_dict_affinity))

boltz_pred_path = definitions.score_path / "boltz_pred_probabilities.json"
boltz_aff_path = definitions.score_path / "boltz_pred_affinities.json"
if boltz_aff_path.exists():
    with open(boltz_pred_path, "r") as f:
        current_data = json.load(f)
    complete_pred_data, changes = update_dict(
        current_data, benchmark_data_dict_probability
    )
    print(f"{changes} changes to the pred dict")
else:
    complete_pred_data = benchmark_data_dict_probability
with open(boltz_pred_path, "w") as fp:
    json.dump(complete_pred_data, fp)


if boltz_aff_path.exists():
    with open(boltz_aff_path, "r") as f:
        current_data = json.load(f)
    complete_aff_data, changes = update_dict(current_data, benchmark_data_dict_affinity)
    print(f"{changes} changes to the aff dict")
else:
    complete_aff_data = benchmark_data_dict_affinity
with open(boltz_aff_path, "w") as fp:
    json.dump(complete_aff_data, fp)

print(len(complete_aff_data))
print(len(complete_pred_data))
