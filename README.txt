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

WEEK 1 - WEEK 2 

The first part of the project is analysing raw RNA-seq data.

PIPELINE:

# It will require use of a SLURM script submission.
# You will also need the Daphnia magna reference genome and related genome annotation from NIES.
# You will be required to have a metadata file.

1) Download SRR files from NCBI via shell script.
2) Quality control using fastqc via command-line, and examine using multiqc.
3) Trimming adapters and poor-quality bases from sequence start.
4) Using Conda, running RSEM genome preparation with STAR, then aligning the data to the NIES's Daphnia magna reference genome.
5) Download counts, merge them into a matrix, and run DESeq script in R.
6) Get significant (p=0.05) differentially expressed genes for each comparison. 

Then, I performed a reciprocal BLAST search for DNMT1, DNMT3a and TET to find the gene names.