#!/bin/bash

#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20 
#SBATCH --mem=200G 
#SBATCH --time=1-02:30:00
#SBATCH --output=joblogs/job.indexGenome.%J.log 
#SBATCH --error=joblogs/job.indexGenome.%J.log
#SBATCH -J "indexGenome" 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=corrinne@iastate.edu 

target=AD1_AtDt.transcripts.fa

./kallisto index -i ${target%%.*}.idx $target
