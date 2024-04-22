#!/bin/bash

#SBATCH --job-name=array
#SBATCH --array=1-56
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20 
#SBATCH --mem=200G 
#SBATCH --time=1-02:30:00
#SBATCH --output=joblogs/job.cotton.%A_%a.log
#SBATCH --error=joblogs/job.cotton.%A_%a.log
#SBATCH -J "cotton-array" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 

rDir=/work/LAS/jfw-lab/corrinne/TM1/trimmed_reads
genome=/work/LAS/jfw-lab/corrinne/TM1/AD1_AtDt.idx
thr=20

file1=$(ls -1 $rDir/*.R1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls -1 $rDir/*.R2.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

name=$(basename $file1 .R1.fq.gz)

mkdir -p mapping
cd /work/LAS/jfw-lab/corrinne/TM1/mapping

./../kallisto quant -t $thr -i $genome -o $name $file1 $file2
