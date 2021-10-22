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
tmp_dir=$(mktemp -d)

echo ${tmp_dir}

echo $protein $tmp_dir
cp ${protein} ${tmp_dir}
echo active_site_ligand $tmp_dir
cp ${active_site_ligand} ${tmp_dir}
echo ligands_to_dock $tmp_dir
cp ${ligands_to_dock} ${tmp_dir}


protein_name=$(basename ${protein})
active_site_name=$(basename ${active_site_ligand})
ligand_name=$(basename ${ligands_to_dock})

echo ${protein_name}
echo $active_site_name
echo $ligand_name
/scratch/gutermuth/software/JamdaPreprocessing_0.9.1/JamdaPreprocessing -p ${tmp_dir}/${protein_name} -l ${tmp_dir}/${active_site_name} -o ${tmp_dir}/${protein_name}_prep.pdb

/scratch/gutermuth/software/JamdaDocker_0.9.0/JamdaDocker initialize -p ${tmp_dir}/${protein_name}_prep.pdb -l ${tmp_dir}/${active_site_name} -o ${tmp_dir}/${protein_name}.receptor

/scratch/gutermuth/software/JamdaDocker_0.9.0/JamdaDocker dock -r ${tmp_dir}/${protein_name}.receptor -m ${tmp_dir}/${ligand_name} --proteinOutput 0 -o ${outputdir} -v 3 --rmsd_to ${active_site_ligand}

cp ${tmp_dir}/${protein_name}_prep.pdb ${outputdir}