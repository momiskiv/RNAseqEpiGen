  ______ _             _  ______          _           _             ______            _           _                                             _____      _                       _   _          
|  ___(_)           | | | ___ \        (_)         | |            |  _  \          | |         (_)                                           |  ___|    (_)                     | | (_)         
| |_   _ _ __   __ _| | | |_/ / __ ___  _  ___  ___| |_   ______  | | | |__ _ _ __ | |__  _ __  _  __ _   _ __ ___   __ _  __ _ _ __   __ _  | |__ _ __  _  __ _  ___ _ __   ___| |_ _  ___ ___ 
|  _| | | '_ \ / _` | | |  __/ '__/ _ \| |/ _ \/ __| __| |______| | | | / _` | '_ \| '_ \| '_ \| |/ _` | | '_ ` _ \ / _` |/ _` | '_ \ / _` | |  __| '_ \| |/ _` |/ _ \ '_ \ / _ \ __| |/ __/ __|
| |   | | | | | (_| | | | |  | | | (_) | |  __/ (__| |_           | |/ / (_| | |_) | | | | | | | | (_| | | | | | | | (_| | (_| | | | | (_| | | |__| |_) | | (_| |  __/ | | |  __/ |_| | (__\__ \
\_|   |_|_| |_|\__,_|_| \_|  |_|  \___/| |\___|\___|\__|          |___/ \__,_| .__/|_| |_|_| |_|_|\__,_| |_| |_| |_|\__,_|\__, |_| |_|\__,_| \____/ .__/|_|\__, |\___|_| |_|\___|\__|_|\___|___/
                                      _/ |                                   | |                                           __/ |                  | |       __/ |                               
                                     |__/                                    |_|                                          |___/                   |_|      |___/
  __       __     
 |_ _|_ _  / _|___ 
  | || ' \|  _/ _ \
 |___|_||_|_| \___/
                   
AUTHOR: Naomi Musto
TITLE: MBiolSci Biological Sciences (Genetics) specialised in Bioinformatics
UNIVERSITY: University of Leicester, University Rd, Leicester, LE1 7RH
DEPARTMENT: Department of Genetics, Genomics and Cancer Sciences
EMAIL: nm471@student.le.ac.uk
SECOND EMAIL: naomimusto.skiv@gmail.com
SUPERVISOR: Dr Hollie Marshall
  ___     _           
 |_ _|_ _| |_ _ _ ___ 
  | || ' \  _| '_/ _ \
 |___|_||_\__|_| \___/
                      

Hello there! I'm Naomi Musto, a MBiolSci Genetics student at the University of Leicester.
I have created this Git repo to showcase the progress of my final project.

Our main research question is whether ancestrale environemnts of specific temporal population affect the expression of DNA methylatation-related genes.
Epigenetic mechanisms act as a means ot regulate gene expression in response to environmental cues. Ancestral exposure to environemntal stress could prime the offspring to mount a beneficial epigenetic response to the stressor.
My aim is to identify whether the expression of DNA methylation genes differs between temporal population affected by disparate environmental conditions.

The data I was given entails raw RNA-seq data from 16 clonal lines, from three temporal populations (including their biological and technical replicates): recovery, pesticide, eutrophic.

The project's main aims are:

 * Analysis of RNA-seq data to assess expression of DNMT1 and DNMT3a homologs in different temporal populations of D. magna. 
 * Reciprocal BLAST to identify TET, another gene involved in DN methylation.
 * Gene Ontology analysis of the differentially expressed genes between the temporal populations.
 * Machine learning subproject using the methylation data to predict phenotype.

  ___ _  _   _      ___           
 | _ \ \| | /_\ ___/ __| ___ __ _ 
 |   / .` |/ _ \___\__ \/ -_) _` |
 |_|_\_|\_/_/ \_\  |___/\___\__, |
                               |_|

WEEK 1 - WEEK 2 - WEEK 3

The first part of the project is analysing raw RNA-seq data (See RNA_seq directory).

PIPELINE:

# It will require use of a SLURM script submission.
# You will also need the Daphnia magna reference genome and related genome annotation from NIES.
# You will be required to have a metadata file.

1) Download SRR files from NCBI via shell script.
2) Quality control using fastqc via command-line, and examine using multiqc.
3) Trimming adapters and poor-quality bases from sequence start.
4) Using Conda, running RSEM genome preparation with STAR, then aligning the data to the NIES's Daphnia magna reference genome.
5) Download counts, merge them into a matrix, and run DESeq script in R.
6) Get significant (p=0.05) differentially expressed genes for each comparison (PvE, RvE, RvP) -- see DE_Analysis directory.

Then, I performed a reciprocal BLAST search for DNMT1, DNMT3a and TET to find the gene names (see Reciprocal_BLAST directory).
Main results:
- DNMT1 = Daphnia magna DNA (cytosine-5)-methyltransferase 1, (LOC116925494) from variants X1-3
- DNMT3a = Daphnia magna DNA (cytosine-5)-methyltransferase 3B (LOC116922472 and LOC116923817) from variants X1,2 and 8
- TET = Daphnia magna DNA N6-methyl adenine demethylase (LOC116917382) from variants X1-3

I was unsure about the results I got from DNMT3a, but it appears that it is highly similar to DNMT3b.
When trying to plot the identified genes, there were not present in any of the differentially expressed genes in the three comparisons. DNMT1 was barely expressed, with just a handful of samples having around 2-3 counts of the gene. DNMT3b variants X1-2 had no expression at all, while variant X8 had relatively higher expression than the previous (max of 38 counts). TET, on the other hand, had mostly low expression (3-11 counts), with the highest count being 43. However, in respect to the rest of the expressed genes within the analysis, these counts are incredibly small and not statistically powerful to extract any conclusions. 
On a more positive note, there seems to be an equal amount of isoforms of DNMTs as found in the literature(Chaturverdi et al., 2023), despite using different genome assemblies (they used LRV0_1, while I am using NIES).

   ___  ___      ___          _    _                  _       _             _         _    
  / __|/ _ \ ___| __|_ _  _ _(_)__| |_  _ __  ___ _ _| |_    /_\  _ _  __ _| |_  _ __(_)___
 | (_ | (_) |___| _|| ' \| '_| / _| ' \| '  \/ -_) ' \  _|  / _ \| ' \/ _` | | || (_-< (_-<
  \___|\___/    |___|_||_|_| |_\__|_||_|_|_|_\___|_||_\__| /_/ \_\_||_\__,_|_|\_, /__/_/__/
                                                                              |__/         

WEEK 4

After spending some time exploring the relevant literature surrounding D. magna and temporal population studies, I attempted a GO enrichment analysis using the differentially expressed gene data from my RNA-seq final analysis in combination with the gene ontology data released by NIES.

PIPELINE (see GO_enrichment_analysis):
1) Created a dataframe from the gene ontology file from NIES containing GID, GO, EVIDENCE, as well as as another dataframe with the relevant gene information (GID and SYMBOL).
2) Created a custom D. magna database using the makeOrgPackage function.
3) Extracted gene names from the objects containing significantly differentially expressed genes across the comparisons. 
4) Performed an enrichGO for biological processes (PvE and RvE = p-value cutoff: 0.01, q-value cutoff: 0.05; RvP = p-value cutoff: 0.2, q-value cutoff: 1), plotted the results and saved dataframe as csv for later use. 

RvP turned to be a quite weak comparison, which aligns with the RNA-seq data. As for the other two, they provided quite interesting results!
- PvE: Most differentially expressed genes were correlated ot translation, cytoplasmic translation and ribosomal biogenesis.
- RvE: Most differentially expressed genes were also related to translation, cytoplasmic translation and ribosomal biogenesis, but there were most genes involved in lipi transport and localisation.

ANALYSIS IS FINALLY COMPLETE! I will focus on completing the writing for chapter one, and I will update shortly once I have made a start on the second part of my project.
