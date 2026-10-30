# BALLROOM A/S generation, docking, and analysis

This repository supports the generation of the BALLROOM A and BALLROOM S benchmarks, docking workflows, and analysis of benchmark and docking results. The workflows use Python scripts alongside institutional databases, external software, and data configured for the environment where they are run.

## Repository layout

Code is organized under `src/`:

- `src/generation/` — benchmark generation, filtering, MMP processing, and contradiction analyses.
- `src/docking/` — preparation of docking inputs, cluster scripts, score extraction, and evaluation.
- `src/analysis/` — publication-related plots and analyses.
- `src/utils/` — shared definitions, database access, and helper functions.
- `tests/` — unit and integration tests.

## Requirements and setup

- Python 3.13 or newer and [uv](https://docs.astral.sh/uv/).
- Access to the ChEMBL and ActivityDB PostgreSQL databases on a given host for database-backed generation and tests.
- Depending on the workflow, external data and software such as Siena, UNICON, AutoDockFR, JAMDA, and docking software. These are managed separately from the Python environment. ZBH software (including Siena, UNICON, and JAMDA) is available through the [SoftwareServer](software.zbh.uni-hamburg.de). Cluster workflows use an SGE environment and its scheduler.

From the repository root, create the local configuration file from the template:

```bash
cp src/ballroom.env.template src/ballroom.env
```

Edit `src/ballroom.env` with valid database credentials and local paths. The template lists the expected names; replace every `change` value with a real value. Do not commit credentials. The configuration is loaded by `src/utils/definitions.py` and includes paths for external tools, benchmark data, the PDB mirror, and docking outputs.

Run Python commands from the repository root. `uv run` creates/updates the project environment and installs the Python dependencies declared in `pyproject.toml`. For example:

```bash
uv run pytest
```

## Workflows

### Generate benchmarks

[`src/generation/create_ballroomas_complete.py`](src/generation/create_ballroomas_complete.py) is the end-to-end generation script. It builds the preliminary benchmark, derives the crystal-based and MMP-based pair benchmarks, downloads/prepares structures, and calculates Siena clusters and summary statistics. It requires working database access, configured external tools, and the relevant input data.

The script can be invoked from the repository root with:

```bash
PYTHONPATH=. uv run python src/generation/create_ballroomas_complete.py
```

**Output note:** each run removes and recreates `gutermuth/benchmark_directory/benchmark_name` in the repository. Preserve any existing results there before running the script again.

Other generation scripts provide individual steps but are not necessary to evoke, including:

- [`make_benchmark_final.py`](src/generation/make_benchmark_final.py) — build the preliminary activity benchmark.
- [`calc_activity_pairs.py`](src/generation/calc_activity_pairs.py) and [`create_mmp_directory.py`](src/generation/create_mmp_directory.py) — create crystal-based and MMP-based activity pairs.
- [`download_benchmark.py`](src/generation/download_benchmark.py) and [`finalize_siena_benchmark.py`](src/generation/finalize_siena_benchmark.py) — prepare structure and pocket data.
- [`check_pair_benchmark_contradictions.py`](src/generation/check_pair_benchmark_contradictions.py) and [`contradictions_complete_chembl.py`](src/generation/contradictions_complete_chembl.py) — contradiction analyses.

### Prepare and run docking

[`src/docking/create_all_dockings.py`](src/docking/create_all_dockings.py) coordinates the local preparation steps: selecting receptor/ligand combinations, preparing files, and writing docking input tables. It expects benchmark files in the paths configured by `base_benchmark_path`, existing score files, a PDB mirror, and the required preparation tools. Its outputs are under `docking_path/general_dockings`.

After preparation, the docking jobs are run separately. The repository includes scripts for DOCK, Vina, and JAMDA; jobs are submitted on the SGE cluster with `qsub`. The shell scripts are configured for the cluster environment and include site-specific paths. When adapting them, update executable, input, output, working-directory, and email settings. For scripts that read one docking-table row per task, set the scheduler task range to the number of rows in that table. Boltz runs use a separate repository.

Use the `extract_*` scripts in `src/docking/` to collect method outputs. [`evaluate_benchmark_dockings.py`](src/docking/evaluate_benchmark_dockings.py) compares methods on both benchmarks; it expects the benchmark CSVs and score files (including DOCK, Vina, JAMDA, Boltz, and neural-network outputs) in the configured locations.

### Analysis

Scripts in [`src/analysis/`](src/analysis/) produce publication-related results and plots. For example, `evaluate_rmsd_only.py` evaluates RMSD-filtered performance, `make_own_performance_plots.py` creates performance plots, and `molecule_plots.py` creates molecule figures. These scripts generally expect benchmark and/or score outputs to have already been generated; see each script for its specific inputs and output paths.

## Configuration and output paths

Paths are defined centrally in [`src/utils/definitions.py`](src/utils/definitions.py). In particular:

- `base_benchmark_path` is expected to contain `presubmission_ballroom/ballrooma` and `presubmission_ballroom/ballrooms`, each with a `complete_data.csv` for docking evaluation.
- `docking_path` holds `general_dockings/` and `scores/`.
- The shared analysis output path is `gutermuth/output/` under the repository.

Output roots vary by workflow; some generation scripts define paths directly in the script. Review the relevant paths when setting up a run, particularly if you want to retain existing results.

## Tests

Run the test suite from the repository root:

```bash
uv run pytest
```

The test suite includes database integration tests that query both databases and a generation integration test that uses configured paths and external data/tools. The generation test recreates `base_benchmark_path/automatic_test`; use a location where replacing that test output is acceptable.


# Resulting datasets
The resulting datasets can be found on the [FDR](https://www.fdr.uni-hamburg.de/record/21498).