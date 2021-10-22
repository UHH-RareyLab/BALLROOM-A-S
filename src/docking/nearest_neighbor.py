from pathlib import Path

import pandas as pd

import src.utils.definitions as definitions

python_executable = Path(
    "/scratch/gutermuth/anaconda3/envs/activityfinderanalysis/bin/python"
)
model_file = Path(
    "/work/stud2019/tharren/share/BallroomMLModel/NaomiML/tools/BallroomAnalysis/models.py"
)

benchmarks = definitions.benchmarks
smiles_dict = {}

for benchmark in benchmarks:
    complete_data = pd.read_csv(benchmark / "complete_data.csv")
    print(complete_data.shape)
    for index, row in complete_data.iterrows():
        target = row["Target"]
        if target not in smiles_dict.keys():
            smiles_dict[target] = {}
        id1 = row["ID1"]
        id2 = row["ID2"]
        smiles1 = row["Smiles1"]
        smiles2 = row["Smiles2"]
        if not id1 in smiles_dict[target]:
            smiles_dict[target][id1] = smiles1
        if not id2 in smiles_dict[target]:
            smiles_dict[target][id2] = smiles2

smiles_data = []
for target in smiles_dict:
    for ligand in smiles_dict[target]:
        smiles_data.append([ligand, smiles_dict[target][ligand], target])
data = pd.DataFrame(
    smiles_data, columns=["Ligand", "canonical_smiles", "target.chembl_id"]
)

output_path = definitions.docking_output_path / "nearest_neighbor.csv"
data.to_csv(output_path, index=False)
