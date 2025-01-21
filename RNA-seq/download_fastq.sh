#!/bin/bash

# Customise for use: mail, account, time

#SBATCH --job-name=downloading_data
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nm471@student.le.ac.uk
#SBATCH --account=evo-epi

### Include when testing; testing passed 14/01/25
# SBATCH --partition=devel

# Final output dir
# Can customise name but change -O in the loop
# Delete if you have another dir of choice, change -O in loop
mkdir -p fastq_mine

# SRA path from config
# Where your SSR files will be; customise
SRA_CACHE="/scratch/evo-epi/nm471/Project/cache/sra"

### Run script in the working directory it was submitted in
cd $SLURM_SUBMIT_DIR

# Load modules; updated on 14/01/25
module load gcc/12.3.0-yxgv2bl
module load sratoolkit/3.0.0-5fetwpi

# Download data from NCBI
# txt needs to be in working dir
prefetch --option-file Accessions.txt 

# Make SRR files into fastq
for file in ${SRA_CACHE}/SRR* 
do
    base_name=$(basename "$file") # need just basename or crashes
    fasterq-dump --split-files "$file"  -O fastq_mine 
done
