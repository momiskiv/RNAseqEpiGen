#!/bin/bash

# Customise for use: mail, account, time

#SBATCH --job-name=alignment
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=14:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nm471@student.le.ac.uk
#SBATCH --account=evo-epi

### Include when testing
# SBATCH --partition=devel

### Run script in the working directory it was submitted in
# Make sure it has: trimmed .fq.gz files, .fa file and .gff file
cd $SLURM_SUBMIT_DIR

### Make a conda environment BEFORE RUNNING
#conda activate
#conda update -n base -c defaults conda

#conda create -n alignment python=3.9 -y
#conda activate alignment

#conda install -n alignment conda conda-libmamba-solver -y
#conda config --set solver libmamba

#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge

#conda install star samtools htseq -y

# Customise
source /scratch/evo-epi/nm471/miniconda3/bin/activate alignment

### Genome prep 
# D. magna genome = 161.5 Mb 
STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles ./NIES_genome.fna \
--sjdbGTFfile ./NIES_genomic.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases 12

### Alignment
for file in $(ls *_trim_1.fq.gz)
do
    base=$(basename $file "_trim_1.fq.gz")
    STAR \
    --runThreadN 24 \
    --genomeDir ./ \
    --readFilesCommand gunzip -c \
    --readFilesIn ${base}_trim_1.fq.gz ${base}_trim_2.fq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${base}
done


### Counting reads
for file in $(ls *Aligned.sortedByCoord.out.bam)
do
    base=$(basename $file "Aligned.sortedByCoord.out.bam")
    htseq-count --format=bam  \
    ${base}Aligned.sortedByCoord.out.bam \
    ./NIES_genomic.gtf \
    > ${base}.counts
done