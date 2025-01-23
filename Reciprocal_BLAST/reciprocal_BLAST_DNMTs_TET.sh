#!/bin/bash

# Customise for use: mail, account, time

#SBATCH --job-name=reciprocal_BLAST_DNMTs_TET
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=5:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=nm471@student.le.ac.uk
#SBATCH --account=evo-epi

#-------------------------------------------------------------------
# DNMT blast
#-------------------------------------------------------------------

# Generate conda env 
# conda create -n agat_env
# conda activate agat_env
# conda install -c bioconda agat

# Modules
module load gcc/12.3.0-yxgv2bl
module load blast-plus/2.13.0-5o3kbvq

# Get protein seqs of all Dmagna genes
agat_sp_extract_sequences.pl --gff NIES_genomic.gtf -f NIES_genome.fna -p -o dmagna_proteins.fa

makeblastdb -in dmagna_proteins.fa -parse_seqids -dbtype prot

# Navigate to http://v2.insect-genome.com/Pcg
# In advanced search, specify gene acronym, then bulk download all PROTEIN .fasta
# for DNMT1 and DNMT3a for 321 and 110 species respectively 
# for TET >1000 sequences from >700 species

#BLAST variables
EVALUE=1e-3
MAX_TARGETS=5

# BLAST database creation and reciprocal searches
for gene in dnmt1 dnmt3a tet; do
    [ -f ${gene}.fa ] || { echo "Error: ${gene}.fa not found"; exit 1; } #check file existance
    
    echo "Processing ${gene}..."
    
    makeblastdb -in ${gene}.fa -parse_seqids -dbtype prot

done

# Parallel BLAST
echo "Parallel BLAST searches..."
parallel -j $SLURM_CPUS_PER_TASK blastp -query dmagna_proteins.fa -db {}.fa -evalue $EVALUE -max_target_seqs $MAX_TARGETS -outfmt 6 -out dmagna_{}.txt ::: dnmt1 dnmt3a tet

# Parallel reciprocal BLAST
echo "Reciprocal BLAST searches..."
parallel -j $SLURM_CPUS_PER_TASK blastp -query {}.fa -db dmagna_proteins.fa -evalue $EVALUE -max_target_seqs $MAX_TARGETS -outfmt 6 -out dmagna_{}_reciprocal.txt ::: dnmt1 dnmt3a tet

echo "All BLAST searches completed."

# Cleanup intermediates
rm *.fa.pdb *.fa.pjs *.fa.pog *.fa.pos *.fa.ptf *fa.pto *.fna.index.dir *.fna.infex.pag *.agat.log





