#! /bin/bash
#$ -j y
#$ -S /bin/bash
#$ -q hpc.q,8c.q,16c.q,40c.q,64c.q,32c.q
#$ -m a -M gutermuth@zbh.uni-hamburg.de
#$ -N BallroomDock
#$ -o /scratch/gutermuth/all_benchmark_dockings/presubmission_ballroom/cluster_output
#$ -e /scratch/gutermuth/all_benchmark_dockings/presubmission_ballroom/cluster_output
#$ -wd /scratch/gutermuth/all_benchmark_dockings
ulimit -c 0
### This is a cluster script I got from Christiane/Flachsi for cluster submission of a jamda docking

echo ${SGE_TASK_ID}

docking_file=/scratch/gutermuth/all_benchmark_dockings/presubmission_ballroom/general_dockings/docking_df_dock.csv
data=$(head -n ${SGE_TASK_ID} ${docking_file} | tail -n 1)
echo $data

/scratch/gutermuth/anaconda3/envs/bigboy/bin/python /work/gutermuth/bob_analyse_skripte/src/docking/execute_DOCK_dockings.py $data

