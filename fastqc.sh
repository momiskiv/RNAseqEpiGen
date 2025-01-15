#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=05:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nm471@student.le.ac.uk
#SBATCH --account=evo-epi

### Include when testing
# SBATCH --partition=devel

### Run script in the working directory it was submitted in
cd $SLURM_SUBMIT_DIR

# Load modules; updated 15/01/25
module load gcc/12.3.0-yxgv2bl
module load py-multiqc/1.14-ateq76h
module load fastqc/0.12.1-hkgpcde

# Run fastqc to look at the quality of the data
for file in $(ls *.fastq)
do
	fastqc -t 24 ${file}
done

# Run multiqc to make a nice report
multiqc ./