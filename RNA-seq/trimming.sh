#!/bin/bash

# Customise for use: mail, account, time

#SBATCH --job-name=trimming
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=03:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nm471@student.le.ac.uk
#SBATCH --account=evo-epi

### Include when testing
# SBATCH --partition=devel

# Directory for trimmed; otherwise report with both
mkdir -p trimfq

### Run script in the working directory it was submitted in
cd $SLURM_SUBMIT_DIR

# Load modules; updated 15/01/2025
module load gcc/12.3.0-yxgv2bl
module load py-multiqc/1.14-ateq76h
module load fastqc/0.12.1-hkgpcde
module load trimmomatic/0.39-tfj436w

# Trim the reads based on the fastqc file
# Compressing trimmed files for storage
for file in $(ls *_1.fastq)
do
	base=$(basename ${file} "_1.fastq")

	trimmomatic PE -threads 24 ${base}_1.fastq ${base}_2.fastq \
    trimfq/${base}_trim_1.fq.gz trimfq/${base}_unpaired_1.fq.gz \
	trimfq/${base}_trim_2.fq.gz trimfq/${base}_unpaired_2.fq.gz \
    ILLUMINACLIP:illumina_adapters_new.fa:2:30:10 MINLEN:50 HEADCROP:10
done 

# Run fastqc again to see what has changed
cd trimfq

for file in $(ls *trim*.fq.gz)
do
	fastqc -t 24 ${file}
done

# Run multiqc to make a nice report
multiqc ./
