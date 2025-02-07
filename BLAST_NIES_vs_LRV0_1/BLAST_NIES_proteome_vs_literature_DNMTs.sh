#!/bin/bash

NIES_PROTEOME="/scratch/evo-epi/nm471/NIES_Dmagna_Genome/dmagna_proteins.fa"
DNMT1="/scratch/evo-epi/nm471/NIES_Dmagna_Genome/DNMT1.fasta" 
DNMT3a="/scratch/evo-epi/nm471/NIES_Dmagna_Genome/DNMT3A.fasta"
OUTPUT_DIR="/scratch/evo-epi/nm471/NIES_Dmagna_Genome"

# short script to retrieved and BLAST sequences, saved as txt in selected output directory

# Modules; make sure it is updated to current version before use
module load gcc/12.3.0-yxgv2bl
module load blast-plus/2.13.0-5o3kbvq

echo "RETRIEVING SEQUENCES FROM NCBI...."

# Change protein to mine for specific use

efetch -db protein -id "KAK4015743.1,KAK4015744.1,KAK4015745.1" -format fasta > $DNMT1
efetch -db protein -id "KAK4005979.1" -format fasta > $DNMT3a

echo "ALL SEQUENCES EXTRACTED....."
echo "CREATING DB..."

# Can customise DB name if needed

makeblastdb -in $NIES_PROTEOME -parse_seqids -dbtype prot -out $OUTPUT_DIR/nies_proteome_db

EVALUE=0.05

echo "STARTING BLAST......."

# Edit DB name if changed previously
# outfmt 6 = tabular 

blastp -query $DNMT1 -db $OUTPUT_DIR/nies_proteome_db -evalue $EVALUE -outfmt 6 -out $OUTPUT_DIR/blastp_dnmt1_results.txt 
blastp -query $DNMT3a -db $OUTPUT_DIR/nies_proteome_db -evalue $EVALUE -outfmt 6 -out $OUTPUT_DIR/blastp_dnmt3a_results.txt 

echo "COMPLETED BLAST.....!"
