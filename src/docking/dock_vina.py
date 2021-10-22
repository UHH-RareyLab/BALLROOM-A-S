import argparse
import json
from pathlib import Path

from vina import Vina

parser = argparse.ArgumentParser()
parser.add_argument("protein", type=str, help="path to the protein")
parser.add_argument("ligand", type=str, help="path to the ligand")
parser.add_argument("box", type=str, help="path to the box file")
parser.add_argument("output", type=str, help="output directory to write prepared")

args = parser.parse_args()

seed = 42
# start_path = Path("/scratch/gutermuth/all_benchmark_dockings/testing/CHEMBL1075323/2GQP_vina_prep")
start_path = Path(
    "/scratch/gutermuth/all_benchmark_dockings/testing/4LOH_1SY_A_401/4LOH_vina_prep"
)

with open(args.box) as infile:
    d = json.load(infile)
    center = d["CENTER"]
    box_size = d["BOX_SIZE"]

box_size[0] = max(22, box_size[0] + 10)
box_size[1] = max(22, box_size[1] + 10)
box_size[2] = max(22, box_size[2] + 10)

print("Using box size:", box_size)
print("Using random seed:", seed)
# receptor = [x for x in start_path.glob("*.pdbqt") if len(x.stem) == 4][0]
# ligand = [x for x in start_path.glob("*.pdbqt") if len(x.stem) != 4][0]
receptor = Path(args.protein)
ligand = Path(args.ligand)
output = Path(args.output)
print(f"receptor to use: {receptor}")
print(f"ligand to dock: {ligand}")
print(f"Seed used: {seed}")
print(f"Output directory: {output}")
if not output.exists():
    output.mkdir()
results = output / (ligand.stem + ".pdbqt")
v = Vina(sf_name="vina", cpu=1, seed=seed)
v.set_receptor(receptor.as_posix())
v.set_ligand_from_file(ligand.as_posix())
v.compute_vina_maps(center=center, box_size=box_size)
v.dock(exhaustiveness=32, n_poses=32)
v.write_poses(results.as_posix(), n_poses=32, energy_range=float("inf"), overwrite=True)
