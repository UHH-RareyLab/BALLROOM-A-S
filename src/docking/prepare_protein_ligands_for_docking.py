import os
import subprocess

import src.utils.definitions as definitions


def check_returns(process, input):
    """
    Helperfunction that handles subprocess processes and saves their output to the
        loggingfile or raises Error
    :param process: subprocess process that has been executed
    """
    if process.returncode != 0:
        print("Wrong going stuff")
        print(input)
        print(process.stdout.decode())
        print(process.stderr.decode())
        print("end wrong stuff")


def prepare_benchmarks(rerun=False):
    for target in definitions.docking_output_path.glob("*"):
        if target.is_file():
            continue
        print(f"Processing target {target}")
        pdbs = [x for x in target.glob("*.pdb")]
        for unprocessed_pdb in pdbs:
            possible_smiles = target / "to_dock.smi"
            if possible_smiles.exists():
                input_ligands = possible_smiles
            else:
                raise RuntimeError
            crystal_ligands = [
                x
                for x in target.glob("*.sdf")
                if x.name.startswith(unprocessed_pdb.stem)
            ]
            assert len(crystal_ligands) == 1
            crystal_ligand = crystal_ligands[0]
            vina_output_dir = target / (unprocessed_pdb.stem + "_vina_prep")
            processed_dir = target / (unprocessed_pdb.stem + "_prepared")
            prepped_ligands_directory = vina_output_dir / "prepped_ligands"
            # Lets define all of the expected output files
            prepared_ligands = input_ligands.parent / (
                input_ligands.stem + "prepared.smi"
            )
            # Now we have set all important variables, now lets get to preparing
            # first, protonate the smiles using unicon
            prepare_ligands = [
                definitions.unicon_path,
                "-i",
                input_ligands,
                "-o",
                prepared_ligands,
                "-t",
                "single",
                "-p",
                "single",
            ]
            if not prepared_ligands.exists() or rerun:
                prepare_ligands_process = subprocess.run(
                    prepare_ligands, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                check_returns(prepare_ligands_process, prepare_ligands)
            with open(prepared_ligands, "r") as f:
                n_molecules = sum(1 for _ in f)
            # Extra directories
            processed_dir.mkdir(exist_ok=True)
            vina_output_dir.mkdir(exist_ok=True)
            prepped_ligands_directory.mkdir(exist_ok=True)

            processed_pdb2 = processed_dir / (unprocessed_pdb.stem + "_prepared.pdb")
            pdbqt_protein = vina_output_dir / (unprocessed_pdb.stem + ".pdbqt")
            # mol2_ligand = vina_output_dir / (input_ligands.stem + ".mol2")
            # pdbqt_ligand = vina_output_dir / (input_ligands.stem + ".pdbqt")
            prepped_ligands_name = prepped_ligands_directory / "prepped_ligands.mol2"
            box_json = vina_output_dir / "box_dimensions.json"
            prepare_protein = [
                definitions.jamda_prep_path,
                "-p",
                unprocessed_pdb.as_posix(),
                "-l",
                crystal_ligand.as_posix(),
                "-o",
                processed_pdb2,
            ]
            convert_receptor_to_pdbqt = [
                (definitions.adfr_suite_path / "bin" / "prepare_receptor").as_posix(),
                "-r",
                processed_pdb2.as_posix(),
                "-o",
                pdbqt_protein.as_posix(),
                "-U",
                "nphs_lps_waters",
            ]
            convert_ligand_to_mol2 = [
                definitions.unicon_path,
                "-i",
                prepared_ligands,
                "-o",
                "prepped_ligands.mol2",
                "-g",
                "3",
                "-s",
                "1",
            ]
            calc_center_of_ligand = [
                definitions.molecule_center_path,
                "--input",
                crystal_ligand,
                "--output",
                box_json,
            ]
            os.chdir(prepped_ligands_directory)
            if not processed_pdb2.exists() or rerun:
                process0 = subprocess.run(
                    prepare_protein, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                check_returns(process0, prepare_protein)
            if not pdbqt_protein.exists() or rerun:
                process1 = subprocess.run(
                    convert_receptor_to_pdbqt,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                check_returns(process1, convert_receptor_to_pdbqt)
            # Dies from qthread so no check here ... seems to work though
            if not (
                n_molecules
                == len([x for x in prepped_ligands_directory.glob("*.mol2")])
                or rerun
            ):
                process2 = subprocess.run(
                    convert_ligand_to_mol2,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                check_returns(process2, convert_ligand_to_mol2)
            if not box_json.exists() or rerun:
                process4 = subprocess.run(
                    calc_center_of_ligand,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                check_returns(process4, calc_center_of_ligand)
            if not (
                n_molecules
                == len([x for x in prepped_ligands_directory.glob("*.pdbqt")])
                or rerun
            ):
                for mol2_ligand in prepped_ligands_directory.glob("*"):
                    convert_mol2_to_pdbqt = [
                        definitions.adfr_suite_path / "bin" / "prepare_ligand",
                        "-l",
                        mol2_ligand.name,
                        "-o",
                        mol2_ligand.stem + ".pdbqt",
                    ]
                    process3 = subprocess.run(
                        convert_mol2_to_pdbqt,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    )
                    check_returns(process3, convert_mol2_to_pdbqt)
