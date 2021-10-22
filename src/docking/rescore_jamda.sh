#! /bin/bash
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -q hpc.q,8c.q,16c.q,40c.q,64c.q,32c.q
ulimit -c 0
### This is a cluster script I got from Christiane/Flachsi for cluster submission of a jamda docking

echo ${SGE_TASK_ID}

docking_file=$1
data=$(head -n ${SGE_TASK_ID} ${docking_file} | tail -n 1)

IFS=, read -r protein active_site_ligand ligands_to_dock outputdir <<< ${data}

echo ${outputdir}
echo $protein
echo $active_site_ligand
echo $ligands_to_dock

/scratch/gutermuth/software/JamdaScorer_1.2.0/JamdaScorer -i $protein -l $active_site_ligand --optimize --protoss --rmsd -o ${outputdir}/opt.sdf
