import json
from pathlib import Path

import pandas as pd

pdbbind_path = Path(
    "/local/gutermuth/backup/bob_analyse_skripte/gutermuth/pdbbind_v2020"
)
ballroom_path = Path(
    "/local/gutermuth/backup/bob_analyse_skripte/gutermuth/benchmark_directory/benchmark_july_25"
)
ballroom_s_path = ballroom_path / "ballrooms"
ballroom_a_path = ballroom_path / "ballrooma"
tobi_path = Path(
    "/local/gutermuth/backup/bob_analyse_skripte/gutermuth/general_data/pdb_ids_tobi_hetcodes.txt"
)
parsed_activity_path = Path(
    "/scratch/gutermuth/master_table/2025_presubmission_activityfinder/parsed_activity_data_raw_enriched.csv"
)


def extract_tobi_file(line):
    pdb = line.split(",")[0]
    hetcode = line.split(",")[1]
    return pdb.upper() + "_" + hetcode


def extract_protein_ligand_name(line):
    pdb = line.split()[0]
    ligand = line.split()[6].replace("(", "").replace(")", "")
    return pdb.upper() + "_" + ligand


def extract_protein_ligands_from_ballroom(benchmark_path):
    data = pd.read_csv(benchmark_path / "complete_data.csv", index_col=0)
    all_ids = set()
    for index, row in data.iterrows():
        id1 = row["ID1"]
        id2 = row["ID2"]
        if not id1.startswith("CHEMBL"):
            all_ids.add(id1)
        if not id2.startswith("CHEMBL"):
            all_ids.add(id2)
    return [str(x) for x in all_ids]


def calc_overlap_between_benchmarks(pdbbind_list, ballroom_list):
    overlapping_complexes = set()
    for complex in pdbbind_list:
        for complex2 in ballroom_list:
            if complex2.startswith(complex):
                overlapping_complexes.add(complex2)
    return overlapping_complexes


tobi_ids = []
with open(tobi_path, "r") as f:
    for line in f.readlines():
        line = line.strip()
        tobi_ids.append(extract_tobi_file(line))
print(tobi_ids)
pl_path = pdbbind_path / "index" / "INDEX_general_PL.2020"
nl_path = pdbbind_path / "index" / "INDEX_general_NL.2020"
pl_ligands = []
nl_ligands = []


with open(pl_path, "r") as f:
    for line in f.readlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        pl_ligands.append(extract_protein_ligand_name(line))

with open(nl_path, "r") as f:
    for line in f.readlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        nl_ligands.append(extract_protein_ligand_name(line))

all_pdbbind = []
all_pdbbind.extend(nl_ligands)
all_pdbbind.extend(pl_ligands)
ballrooms_ligands = extract_protein_ligands_from_ballroom(ballroom_s_path)
ballrooma_ligands = extract_protein_ligands_from_ballroom(ballroom_a_path)
ballrooma_overlap = calc_overlap_between_benchmarks(all_pdbbind, ballrooma_ligands)
ballrooms_overlap = calc_overlap_between_benchmarks(all_pdbbind, ballrooms_ligands)
tobi_overlap = calc_overlap_between_benchmarks(all_pdbbind, tobi_ids)
with open(
    "../identity_to_benchmarks/pdbligands_overlap_pdbbind_ballrooms.json", "w"
) as f:
    json.dump([x for x in ballrooms_overlap], f)
with open(
    "../identity_to_benchmarks/pdbligands_overlap_pdbbind_ballrooma.json", "w"
) as f:
    json.dump([x for x in ballrooma_overlap], f)
with open("../identity_to_benchmarks/pdbligands_overlap_pdbbind_tobi.json", "w") as f:
    json.dump([x for x in tobi_overlap], f)

# parsed_activities = pd.read_csv(parsed_activity_path, index_col=0)
pdb_chembl_dictionary = {}
for chunk in pd.read_csv(parsed_activity_path, index_col=0, dtype=str, chunksize=1000):
    for index, row in chunk.iterrows():
        pdb = row["pdb.pdb_id"]
        chembl_target = row["target.chembl_id"]
        if pdb in pdb_chembl_dictionary.keys():
            if chembl_target not in pdb_chembl_dictionary[pdb]:
                pdb_chembl_dictionary[pdb].append(chembl_target)
        else:
            pdb_chembl_dictionary[pdb] = [chembl_target]


pdbbind_targets = {}
for complex in all_pdbbind:
    pdb_id = complex.split("_")[0]
    if pdb_id in pdb_chembl_dictionary.keys():
        pdbbind_targets[pdb_id] = pdb_chembl_dictionary[pdb_id]
    else:
        pdbbind_targets[pdb_id] = ["Unknown"]

print(pdbbind_targets)
target_count = {}
for pdb in pdbbind_targets:
    targets = pdbbind_targets[pdb]
    for target in targets:
        if target in target_count.keys():
            target_count[target] += 1
        else:
            target_count[target] = 1

sorted_targets = sorted(target_count.items(), key=lambda x: x[1], reverse=True)

for target in sorted_targets:
    print(target)

with open("../identity_to_benchmarks/pdbbind_activityfinder_targets", "w") as f:
    json.dump(pdbbind_targets, f)

with open("../identity_to_benchmarks/pdbbind_chembl_targets_list", "w") as f:
    json.dump(sorted_targets, f)

print(len(ballrooma_overlap))
print(len(ballrooms_overlap))
print(len(ballrooms_ligands))
print(len(ballrooma_ligands))
print(len(ballrooma_overlap) / len(ballrooma_ligands))
print(len(ballrooms_overlap) / len(ballrooms_ligands))
print(len(tobi_overlap) / len(tobi_ids))
print(len(tobi_ids))
print(tobi_overlap)
