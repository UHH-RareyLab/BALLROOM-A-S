from pathlib import Path
import os
from dotenv import load_dotenv
env_path = Path(__file__).resolve().parent.parent / "ballroom.env"
load_dotenv(dotenv_path=env_path)

base_path = Path(__file__).resolve().parent.parent.parent
output_path = base_path / "gutermuth" / "output"
### Fixed definitions
matchings_molecules = [
    "standard inchi-key matches",
    "canonical smiles with chiral information annotated",
    "canonical smiles without chiral information annotated",
    "truncated standard inchi matches if cut after hydrogen connection layer",
    "truncated standard inchi matches if cut after atom connection layer",
]

general_filters = {
    "heavy_atoms_min": 5,
    "heavy_atoms_max": 50,
    "resolution_max": 2.5,
    "dpi_max": 0.42,
    "rfactor_max": 0.45,
}

matchings_proteins = ["Gold", "Silver", "Bronze"]

cutoff_range = 0.2

### Changeable definitions
allowed_activities = ["Ki", "Kd", "IC50", "EC50"]

mm_medium_strict = [
    "standard inchi-key matches",
    "canonical smiles with chiral information annotated",
    "canonical smiles without chiral information annotated",
]
mm_strict = [
    "standard inchi-key matches",
    "canonical smiles with chiral information annotated",
]

### Paths

sienatools_path = Path(os.getenv("sienatools_path"))
siena_path = Path(os.getenv("siena_path"))
unicon_path = Path(os.getenv("unicon_path"))

base_stractable_path = Path(os.getenv("base_stractable_path"))
stractablef_path = base_stractable_path / "StrAcTable_filtered.csv"
stractable_path = base_stractable_path / "StrAcTable.csv"
active_site_data_path = base_stractable_path / "active_site_data.csv"

mirror_base_path = Path(os.getenv("mirror_base_path"))
mirror_path = mirror_base_path / "data" / "structures" / "all" / "pdb"
rdkit_path = Path(os.getenv("rdkit_path"))

### URLS
pdb_download_url = "http://files.rcsb.org/download/"


pair_data_columns = [
    "ID1",
    "ID2",
    "Smiles1",
    "Smiles2",
    "Activity1",
    "Activity2",
    "Difference",
    "Activity_Type",
    "Activity_Unit",
    "Assay",
    "Target",
    "Dockinto1",
    "Dockinto2",
    "Activity_ID_1",
    "Activity_ID_2",
]

# DB shenaniganz
activitydb_username = os.getenv("activitydb_username")
activitydb_password = os.getenv("activitydb_password")
activitydb_databasename = os.getenv("activitydb_name")
chembl_username = os.getenv("chembl_username")
chembl_password = os.getenv("chembl_password")
chembl_databasename = os.getenv("chembl_name")


adfr_suite_path = Path(os.getenv("adfr_suite_path"))
molecule_center_path = Path(os.getenv("molecule_center_path"))
jamda_prep_path = Path(os.getenv("jamda_prep_path"))
base_benchmark_path = Path(os.getenv("base_benchmark_path"))
benchmark_path = base_benchmark_path / "presubmission_ballroom"
ballrooma = benchmark_path / "ballrooma"
ballrooms = benchmark_path / "ballrooms"
benchmarks = [ballrooma, ballrooms]

docking_path = Path(os.getenv("docking_path"))
docking_path.mkdir(parents=True, exist_ok=True)
docking_output_path = docking_path / "general_dockings"
docking_output_path.mkdir(exist_ok=True)

score_path = docking_path / "scores"
score_path.mkdir(exist_ok=True)
