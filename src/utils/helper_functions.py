import json
import shutil

import pandas as pd

from src.generation.check_pair_benchmark_contradictions import check_targets_assays


def update_dict(dict1, dict2):
    changes = 0
    for key in dict2:
        if key in dict1:
            for subkey in dict2[key]:
                if subkey in dict1[key]:
                    if dict1[key][subkey] != dict2[key][subkey]:
                        dict1[key][subkey] = dict2[key][subkey]
                        changes += 1
                else:
                    dict1[key][subkey] = dict2[key][subkey]
                    changes += 1
        else:
            dict1[key] = dict2[key]
    return dict1, changes


def cleanup_pair_benchmark(benchmark_path):
    for target in benchmark_path.glob("*"):
        if target.is_file():
            continue
        pair_data_path = target / "pair_data.csv"
        pair_data = pd.read_csv(pair_data_path)
        if pair_data.shape[0] == 0:
            print(f"removing {target}")
            shutil.rmtree(target)


def extract_pdb_ligand_from_name(name):
    if name.startswith("CHEMBL") or name == "":
        return None
    else:
        return [name.split("_")[0], name]


def get_all_mmp_data_in_benchmark(benchmark_path):
    targets = [x for x in benchmark_path.glob("*") if x.is_dir()]
    mmps = []
    for target in targets:
        for assay in target.glob("*"):
            if assay.is_file():
                continue
            for mmp in assay.glob("*.csv"):
                mmps.append(mmp)
    return mmps


def create_stats_for_benchmark(benchmark_path):
    targets, assays, datapoints, pockets, ligands, pairs = check_targets_assays(
        benchmark_path, []
    )
    ballrooma_stats = {
        "N_Targets": targets,
        "N_assays": assays,
        "N_datapoints": datapoints,
        "N_pockets": pockets,
        "N_ligands": ligands,
        "N_pairs": pairs,
    }
    with open(benchmark_path / "stats.json", "w") as f:
        json.dump(ballrooma_stats, f)


def prop_per_x(x, count):
    """
    Compute the proportion of the counts for each value of x
    """
    df = pd.DataFrame({"x": x, "count": count})
    prop = df["count"] / df.groupby("x")["count"].transform("sum")
    return prop
