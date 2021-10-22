import json
from pathlib import Path

import pandas as pd

latex_dict = {}

replacement_dict = {
    0: "zero",
    1: "one",
    2: "two",
    3: "three",
    4: "four",
    5: "five",
    6: "six",
}


def exchange_int_for_string(number_to_check):
    print(number_to_check)
    if type(number_to_check) != int:
        return number_to_check
    for key, value in replacement_dict.items():
        if key == number_to_check:
            return value
    return number_to_check


def extract_json_stats(path, benchmark="ballrooma"):
    # {"N_Targets": 210, "N_assays": 642, "N_datapoints": 1888, "N_pockets": 850, "N_ligands": 1406, "N_pairs": 1546}
    with open(path, "r") as f:
        data = json.load(f)
    latex_dict[f"{benchmark}_targets"] = data["N_Targets"]
    latex_dict[f"{benchmark}_assays"] = data["N_assays"]
    latex_dict[f"{benchmark}_datapoints"] = data["N_datapoints"]
    latex_dict[f"{benchmark}_pockets"] = data["N_pockets"]
    latex_dict[f"{benchmark}_ligands"] = data["N_ligands"]
    latex_dict[f"{benchmark}_pairs"] = data["N_pairs"]


def generate_latex_definition(
    command, value, percent=False, isstring=False, multiply_by_100=True
):
    if not isstring:
        if not percent:
            value = exchange_int_for_string(value)
            if type(value) == str:
                value = f"{value}" + "\\xspace"
            else:
                value = f"\\num{{{value}}}" + "\\xspace"
        else:
            if multiply_by_100:
                factor = 100
            else:
                factor = 1
            value = f"\\pct{{{factor * value:.2f}}}" + "\\xspace"
    else:
        value = f"{value}\\xspace"
    command = "\\" + str(command)
    command = command.replace("_", "")
    return f"\\newcommand{{{command}}}{{{value}}}\n"


base_path = Path(
    "/work/gutermuth/bob_analyse_skripte/gutermuth/benchmark_directory/presubmission_activityfinder_detailed_contradictions_test"
)
ballrooma_path = base_path / "ballrooma"
ballrooms_path = base_path / "ballrooms"
extract_json_stats(ballrooma_path / "stats.json", "ballrooma")
extract_json_stats(ballrooms_path / "stats.json", "ballrooms")
complete_data_s = pd.read_csv(ballrooms_path / "complete_data.csv")
complete_data_a = pd.read_csv(ballrooma_path / "complete_data.csv")
ballrooma_targets = set([x for x in complete_data_a["Target"]])
ballrooms_targets = set([x for x in complete_data_s["Target"]])
ballrooma_assays = set([x for x in complete_data_a["Assay"]])
ballrooms_assays = set([x for x in complete_data_s["Assay"]])
intersect_targets = ballrooma_targets.intersection(ballrooms_targets)
intersect_assays = ballrooma_assays.intersection(ballrooms_assays)
latex_dict["overlap_targets"] = len(intersect_targets)
latex_dict["overlap_assays"] = len(intersect_assays)

with open("latex_commands.tex", "w") as f:
    for key_value in latex_dict.items():
        f.write(generate_latex_definition(key_value[0], key_value[1]))
