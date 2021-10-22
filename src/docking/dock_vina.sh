#! /bin/bash
#$ -j y
#$ -S /bin/bash
#$ -q hpc.q,8c.q,16c.q,40c.q,64c.q,32c.q
#$ -m a -M gutermuth@zbh.uni-hamburg.de
#$ -o /scratch/gutermuth/all_benchmark_dockings/presubmission_activityfinder_bugfix/cluster_output
#$ -e /scratch/gutermuth/all_benchmark_dockings/presubmission_activityfinder_bugfix/cluster_output
#$ -wd /scratch/gutermuth/all_benchmark_dockings

ulimit -c 0
### This is a cluster script I got from Christiane/Flachsi for cluster submission of a jamda docking

line=$(head -n ${SGE_TASK_ID} /scratch/gutermuth/all_benchmark_dockings/presubmission_ballroom/general_dockings/docking_df_vina.csv | tail -n 1)
echo $line

arr=($line)
echo ${arr[0]}
echo ${arr[1]}
echo ${arr[2]}
echo ${arr[3]}

/scratch/gutermuth/anaconda3/envs/vina_tryaround/bin/python /scratch/gutermuth/all_benchmark_dockings/dock_vina.py ${arr[0]} ${arr[1]} ${arr[2]} ${arr[3]}

echo "Work done"
