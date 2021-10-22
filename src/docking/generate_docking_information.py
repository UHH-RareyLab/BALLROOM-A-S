from pathlib import Path

import pandas as pd

pair_benchmark = Path("/scratch/gutermuth/benchmark_dockings/pair_benchmark")
redocking_data = []
for target in pair_benchmark.glob("*"):
    if target.is_file():
        continue
    ligands = [x for x in target.glob("*.sdf")]
    for ligand in ligands:
        docking_output = target / (ligand.stem + "_docking")
        if not docking_output.exists():
            docking_output.mkdir()
        pdb = target / (ligand.stem.split("_")[0] + ".pdb")
        to_append = [pdb, ligand, ligand, docking_output]
        # print(to_append)
        redocking_data.append(to_append)
df = pd.DataFrame(redocking_data)
# df.to_csv("/scratch/gutermuth/benchmark_dockings/redocking_data.csv", sep=",", index=None, header=None)


pair_benchmark_mmp = Path("/scratch/gutermuth/benchmark_dockings/pair_benchmark_mmp")

mmp_docking_data = []
for target in pair_benchmark_mmp.glob("*"):
    if target.is_file():
        continue
    pair_data_path = target / "pair_data.csv"
    pair_data = pd.read_csv(pair_data_path)
    if pair_data.shape[0] == 0:
        continue
    for assay in target.glob("*_heavy_atoms"):
        for mmp in assay.glob("*.csv"):
            mmp_data = pd.read_csv(mmp)
            smiles = []
            complexes = []
            smiles_file = assay / (mmp.stem + ".smi")
            got_complex_file = False
            for index, row in mmp_data.iterrows():
                name = row["Name"]
                if not name.startswith("CHEMBL"):
                    complexes.append(name)
            for complex in complexes:
                current_ligand_file = target / (complex + "_ligand.sdf")
                current_complex_file = target / (complex.split("_")[0] + ".pdb")
                print(
                    smiles_file.exists(),
                    current_ligand_file.exists(),
                    current_complex_file.exists(),
                )
                print(
                    smiles_file.name,
                    current_ligand_file.name,
                    current_complex_file.name,
                )
                if not current_complex_file.exists():
                    print(current_complex_file)
                    continue
                else:
                    got_complex_file = True
            if not got_complex_file:
                print(assay.name, target.name, complexes, mmp.name)
                raise RuntimeError
