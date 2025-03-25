# GO enrichment analysis 
# Based on tutorial by Sanbomics: https://www.youtube.com/watch?v=JPwdqdo_tRg

setwd("/home/nm471/final_project/github/GO_enrichment_analysis/Differential_Meth_Enrichment")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi")

library(clusterProfiler)
library(AnnotationDbi)
library(readr)
library(dplyr)

#reading gaf file
gaf_data <- read.delim("NIES_gene_ontology.gaf",
                       sep = "\t",
                       comment.char = "!",
                       header = F,
                       col.names = c("!#DB", "GeneID", "Symbol", "Qualifier", "GO_ID", "Reference", "Evidence_Code", "With,From", "Aspect", "Gene_Name", "Gene_Synonym", "Type", "Taxon", "Date", "Assigned_By", "Annot_Ext", "Gene_Product_Form_ID"))

#format for package creation
go_data <- gaf_data %>%
  select(GID = GeneID, GO = GO_ID, EVIDENCE = Evidence_Code) %>%
  distinct()

gene_info <- gaf_data %>%
  select(GID = GeneID, SYMBOL = Symbol) %>%
  distinct()

#I have already created the database before, but do uncomment it and run it if you need to make it (with the changes needed)

# #creating custom D. magna database 
# makeOrgPackage(
#   go = go_data,
#   gene_info = gene_info,
#   version = "0.3",
#   author = "Naomi Musto <nm471@student.le.ac.uk>",
#   maintainer = "Naomi Musto <nm471@student.le.ac.uk>",
#   outputDir = ".",
#   tax_id = "35525",
#   genus = "Daphnia",
#   species = "magna",
#   goTable = "go")

# if need to remake it then do:
# unlink("./org.Dmagna.eg.db", recursive = TRUE, force = TRUE)

#./org.Dmagna.eg.db 

install.packages("/home/nm471/final_project/github/GO_enrichment_analysis/org.Dmagna.eg.db", repos = NULL, type = "source")
library(org.Dmagna.eg.db)

#getting gene names
# 
PvE_diffmeth_genes <- read_csv("/home/nm471/final_project/github/Diff_Meth/PvE_sign_diff_meth_genes.csv", col_select = 1) %>%
  pull(1)

# RvE_diffmeth_genes <- read_csv("/home/nm471/final_project/github/Diff_Meth/RvE_sign_diff_meth_genes.csv", col_select = 1) %>%
#   pull(1)

# RvP_diffmeth_genes <- read_csv("/home/nm471/final_project/github/Diff_Meth/RvP_sign_diff_meth_genes.csv", col_select = 1) %>%
#   pull(1)


PvE_GO_results <- enrichGO(gene = PvE_diffmeth_genes,
                           OrgDb = org.Dmagna.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.5,
                           keyType = "SYMBOL")

as.data.frame(PvE_GO_results)
plot(barplot(PvE_GO_results, showCategory = 10))
write.csv(PvE_GO_results, "PvE_GO_Results.csv")
