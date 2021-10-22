from create_docking_paths import (
    write_docking_dfs,
)
from figure_out_dockings import figure_out_dockings
from prepare_protein_ligands_for_docking import prepare_benchmarks

# figure_out_dockings schreibt den benchmark von local auf scratch und addiert erste docking files
figure_out_dockings()
# prepare protein_ligands_for_docking um die vina_prep sachen zu generieren
prepare_benchmarks()
# create_docking paths um die docking files für dock und vina zu generieren
write_docking_dfs()
# generate docking information auch wieder in die richtung, aber diesmal outputtet das irgendwie nichts
# generate docking inputs ist auch wieder irgendwie dasselbe ohne das es was macht. Was hab ich damals geraucht

# Hier sollte ich nun die bereits existierenden dockings anschauen und nur das delta übrig lassen

# dock jamda mit dem clusterscript  qsub -N BallroomJamda -M gutermuth@zbh.uni-hamburg.de -t 2-1653 -m a /work/gutermuth/bob_analyse_skripte/benchmark_docking_scripts/prepare_and_dock_jamda.sh /scratch/gutermuth/all_benchmark_dockings/presubmission_activityfinder_bugfix/general_dockings/docking_df_jamda.csv
# dock vina um vina zu docken
# generate mw testdate generiert genau das
# dock_dock.sh für die dock dockings

# restart specific dock dockings

# extract boltz
# extract denvis
# extract dock
# extract jamda
# extract_jamda_dockings()
# extract vina
# evaluate_benchmark_dockings um die performance zu berechnen


# extract pdbid ligands altes skript um die ähnlichkeit zu pdbbind zu bestimmen? muss ich mir nochmal anschauen...


# qsub -N BallroomJamda -M gutermuth@zbh.uni-hamburg.de -t 2-1653 -m a /work/gutermuth/bob_analyse_skripte/src/docking/prepare_and_dock_jamda.sh /scratch/gutermuth/all_benchmark_dockings/presubmission_ballroom/general_dockings/docking_df_jamda.csv
