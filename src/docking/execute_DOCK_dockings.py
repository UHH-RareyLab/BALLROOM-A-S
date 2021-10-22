import sys
from pathlib import Path

sys.path.append("/work/gutermuth/dock_scripts")
from cross_docking import CrossDocking
import configparser
import logging
import subprocess
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("folder", type=str, help="path to the receptor folder")

logging.basicConfig(level=logging.DEBUG)
unicon_path = Path("/scratch/gutermuth/software/unicon_1.4.1/unicon")
config_path = Path("/work/gutermuth/dock_scripts/config.ini")
config = configparser.ConfigParser()
config.read(config_path)


def calc_crossdocking(input_path):
    pdbs = [x for x in input_path.glob("*.pdb")]
    assert len(pdbs) == 1
    pdb = [x for x in input_path.glob("*.pdb")][0]
    crystal_ligands = [
        x for x in input_path.glob("*.sdf") if x.name.startswith(pdb.stem)
    ]
    assert len(crystal_ligands) == 1
    crystal_ligand = crystal_ligands[0]
    to_dock = input_path / "to_dockprepared.smi"
    dock_output_path = input_path / "DOCK"
    if not dock_output_path.exists():
        dock_output_path.mkdir()
    prepped_ligand_output = dock_output_path / "to_dockprepared.sdf"
    print(pdb, crystal_ligand, prepped_ligand_output, dock_output_path)
    todo = [str(unicon_path), "-i", to_dock, "-o", prepped_ligand_output, "-g", "3"]
    subprocess.run(todo, check=True)
    cross_docking = CrossDocking(
        str(pdb),
        str(crystal_ligand),
        str(prepped_ligand_output),
        str(dock_output_path),
        config,
    )
    cross_docking.run()


args = parser.parse_args()
folder = args.folder
folder = Path(folder)
real_folder = Path("/work/gutermuth/general_dockings") / folder.name
calc_crossdocking(Path(real_folder))
