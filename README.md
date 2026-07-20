# Ballroom A / S Generation and docking scripts

This repo contains all code that is used for the Ballroom A/S generation, docking and analysis.
The repo is split into four different folders, 
generation for benchmark generation, docking for everything related to docking and scoring,
analysis for analysing both the benchmark and the docking and 
utils for everything that is used by all folders. 

## Python environment
This repo uses uv, you can just uv run from the base repo and everything should install and work.

## Generation 

[create_ballroomas_complete.py](src/generation/create_ballroomas_complete.py) 
Is the main script for generating a new benchmark. 
Most relevant options like the activity difference or linking quality can be directly altered here.

[calc_activity_pairs.py](src/generation/calc_activity_pairs.py) 
Is used to calculate all Molecule/Activity pairs in both benchmarks.

[check_pair_benchmark_contradictions.py](src/generation/check_pair_benchmark_contradictions.py) 
Calculates contradictions between pairs.

[contradictions_complete_chembl.py](src/generation/contradictions_complete_chembl.py) 
Can be used to reproduce the ChEMBL contradiction analysis.

[create_mmp_directory.py](src/generation/create_mmp_directory.py) 
Is used for Ballroom S benchmark generation.

[datasetFilter.py](src/generation/datasetFilter.py) 
Filters each assay dataset.

[download_benchmark.py](src/generation/download_benchmark.py) 
"Downloads" (rather copies) all pdb structures, extracts ligands, etc...

[finalize_siena_benchmark.py](src/generation/finalize_siena_benchmark.py) 
Calculates the siena pocket data.

[generate_timecodes_for_benchmarks.py](src/generation/generate_timecodes_for_benchmarks.py) 
This data can be used to append timecodes to the benchmark if one wants to create a timesplit. Not used in the final publication.

[make_benchmark_final.py](src/generation/make_benchmark_final.py) 
Is used to generate the preliminary benchmark that both BALLROOM A and S are based on.

[make_siena_clusters.py](src/generation/make_siena_clusters.py)
Calculates the connected components in Siena pocket data to find identical/overlapping pockets.

[parse_mmp_output.py](src/generation/parse_mmp_output.py) 
This handles the Output for multiple MMP versions, namely a builtin NAOMI one (not used in the final version) and the MMPDB version.
Also calculates MMP Series. 

[retrieve_bioactivity_info_from_chembl.py](src/generation/retrieve_bioactivity_info_from_chembl.py)
Retrieves additional data from ChEMBL where necessary. 

## Docking

[create_all_dockings.py](src/docking/create_all_dockings.py) is the main script for the docking part. 
This encompasses multiple different things, mainly figuring out which dockings to do ([figure_out_dockings.py](src/docking/figure_out_dockings.py)),
preparing the pdb and ligand files for the dockings ([prepare_protein_ligands_for_docking.py](src/docking/prepare_protein_ligands_for_docking.py))
and writing the different docking csvs so they can be done on the cluster ([create_docking_paths.py](src/docking/create_docking_paths.py)).
Afterwards, you need to manually start all runs for DOCK, Vina, JAMDA and Boltz. 

JAMDA Dockings can be started at the ZBH with the following command:
```
qsub -N BallroomJamda -M user@zbh.uni-hamburg.de -t 1-n -m a /pathto/prepare_and_dock_jamda.sh /pathto/docking_df_jamda.csv
```
JAMDA rescoring can be started at the ZBH with the following command:
```
qsub -N BallroomRescoring -M user@zbh.uni-hamburg.de -t 1-n -m a /pathto/rescore_jamda.sh /pathto/docking_df_jamda.csv
```
DOCK Dockings can be started at the ZBH with the following command:
```
qsub -t 1-n /pathto/dock_dock.sh
```
Vina Dockings can be started at the ZBH with the following command:
```
qsub -t 1-n /pathto/dock_vina.sh
```
Boltz Runs can be started using a different repo.

Afterwards, the extract_ scripts can be used to extract all scores and save them to the score directory to later analyze.

The Nearest Neighbor can be done by generating a docking file using the [nearest_neighbor.py](src/docking/nearest_neighbor.py) file.
Using the paths in the python script a scored output file can be generated from the docking file in the scores subdirectory.

Finally, the [evaluate_benchmark_dockings.py](src/docking/evaluate_benchmark_dockings.py) script can be used to calculate the performance of all methods on both benchmarks.

## Analysis
This part contains all the code used to analyse both benchmarks and all additional analysis in the data.

[analyse_chembl_pair_contradictions.py](src/analysis/analyse_chembl_pair_contradictions.py)
Creates the ChEMBL only pair contradiction analysis for the publication. 

[create_latex_defintions.py](src/analysis/create_latex_defintions.py)
This is my lazy way to generate numbers for the publication. 
Sadly, I cannot generate all numbers that way.

[evaluate_rmsd_only.py](src/analysis/evaluate_rmsd_only.py)
RMSD Performance analysis for the publication.

[make_descriptional_benchmark_plots.py](src/analysis/make_descriptional_benchmark_plots.py)
Additional description plots like activity types for the publication.

[make_own_performance_plots.py](src/analysis/make_own_performance_plots.py)
Make the performance plots for the publication.

[molecule_plots.py](src/analysis/molecule_plots.py)
Create the molecule plots for the publication.

[sunburst.py](src/analysis/sunburst.py)
Create a sunburst plot for each benchmark. 
This is not directly used in the publication, but the numbers calculated are. 

## Utils
Here are some small but important functions and definitions.

[db_interaction.py](src/utils/db_interaction.py)
Handles queries to either the local ChEMBL version or ActivityDB.

[helper_functions.py](src/utils/helper_functions.py)
Some small helper functions.

[definitions.py](src/utils/definitions.py)
Here all relevant definitions, e.g. StrAcTable path, pdb mirror path, tool paths etc. are specifiied. 
All relevant paths, passwords or similar definitions have to be provided using the ballroom.env.template file.
You need to adjust many of these if you want to generate your own benchmark. 
You will also need many tools from the Naomi ChemBioSuite.

## Tests 
Tests can be run using the command uv run pytest