import json

from rdkit import Chem

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


def get_scores_from_mol2_file(mol2_file_path):
    supplier = Chem.ForwardSDMolSupplier(str(mol2_file_path))
    for mol in supplier:
        if mol is None:
            continue
        print(mol.GetNumAtoms())


def print_error_in_log_file(log_file):
    with open(log_file, "r") as f:
        lines = f.readlines()
    print_stuff = False
    for line in lines:
        line = line.strip()
        if line.startswith("raise"):
            print_stuff = True
        if print_stuff:
            print(line)


def get_scores_from_log_file(log_file):
    with open(log_file, "r") as f:
        lines = f.readlines()
    target = None
    scores = {}
    current_molecule = None
    current_score = None
    for line in lines:
        line = line.strip()
        if line.startswith("/work/gutermuth") and target == None:
            target = line.split()[1].split("/")[-1]
            # print(benchmark_type, receptor)
        elif line.startswith("Molecule:"):
            current_molecule = line.split()[1]
            # print(f"current molecule found {current_molecule}")
        elif line.startswith("Grid_Score:"):
            if not current_molecule:
                raise RuntimeError("Found score without molecule")
            if target:
                current_molecule = receptor + "-" + current_molecule
            current_score = float(line.split()[1])
            # print(f"New score of {current_score} found for molecule {current_molecule}")
            if not current_molecule in scores.keys():
                scores[current_molecule] = current_score
            else:
                if current_score < scores[current_molecule]:
                    scores[current_molecule] = current_score
                    # print(f"Better score for molecule {current_molecule} found ")

    return [target, scores]


# logfilepath = Path("/scratch/gutermuth/all_benchmark_dockings/dock_logs/dock_dockings.sh.o451052.100")
logfile_general = definitions.docking_path / "cluster_output"
# logfile_general = Path("/scratch/gutermuth/all_benchmark_dockings/presubmission_activityfinder_bugfix/cluster_output")
logfile_tasknumber = 1349691
logfile_pattern = f"BallroomDock.o{logfile_tasknumber}.*"

# redocking_benchmark_path = Path("/scratch/gutermuth/all_benchmark_dockings/benchmark_dockings/pair_benchmark")
# mmp_benchmark_path = Path("/scratch/gutermuth/mmp_pair_benchmark_dockings")


DOCK_scores = {}
redocking_scores = {}
general_scores = {}
missing_numbers = []

n_errors = 0
n_ges = 0
for logfile in logfile_general.glob(logfile_pattern):
    print(logfile)
    if logfile.name.endswith(".swp"):
        continue
    n_ges += 1
    receptor, scores = get_scores_from_log_file(logfile)
    if len(scores) == 0:
        print(receptor, scores, logfile)
        n_errors += 1
        print_error_in_log_file(logfile)
        missing_number = int(logfile.name.split(".")[-1])
        missing_numbers.append(missing_number)
        continue
    trimmed_receptor = receptor[:-11]
    pdb = trimmed_receptor.split("_")[0]
    ligand = trimmed_receptor[5:]
    complex = pdb + "-" + ligand
    DOCK_scores[complex] = scores
    print(receptor, complex, scores)


print(n_errors, n_ges, n_errors / n_ges)
missing_numbers.sort()

dock_score_path = definitions.score_path / "DOCK.json"

# with open(definitions.score_path / "DOCK.json", "w") as fp:
#    json.dump(DOCK_scores, fp)

if dock_score_path.exists():
    with open(dock_score_path, "r") as f:
        current_jamda_data = json.load(f)
    complete_dock_data, changes = update_dict(current_jamda_data, DOCK_scores)
    print(f"There have been {changes} updates to the dock dict")
else:
    complete_dock_data = DOCK_scores
with open(dock_score_path, "w") as f:
    json.dump(complete_dock_data, f)
