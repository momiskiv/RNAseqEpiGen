# GO enrichment analysis 

setwd("/home/nm471/final_project/github/GO_enrichment_analysis/")

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
gaf_data <- read.delim("GCF_020631705.1_ASM2063170v1.1_gene_ontology.gaf",
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

#creating custom D. magna database 
makeOrgPackage(
  go = go_data,
  gene_info = gene_info,
  version = "0.3",
  author = "Naomi Musto <nm471@student.le.ac.uk>",
  maintainer = "Naomi Musto <nm471@student.le.ac.uk>",
  outputDir = ".",
  tax_id = "35525",
  genus = "Daphnia",
  species = "magna",
  goTable = "go")

# if need to remake it then do:
# unlink("./org.Dmagna.eg.db", recursive = TRUE, force = TRUE)

#./org.Dmagna.eg.db 

install.packages("/home/nm471/final_project/github/GO_enrichment_analysis/org.Dmagna.eg.db", repos = NULL, type = "source")
library(org.Dmagna.eg.db)

#getting gene names

PvE_diffexp_genes <- read_tsv("/home/nm471/final_project/github/RNA-seq/DE_Analysis/sign_diffexp_PvE.tsv", col_select = 1) %>%
  pull(1)

RvE_diffexp_genes <- read_tsv("/home/nm471/final_project/github/RNA-seq/DE_Analysis/sign_diffexp_RvE.tsv", col_select = 1) %>%
  pull(1)

RvP_diffexp_genes <- read_tsv("/home/nm471/final_project/github/RNA-seq/DE_Analysis/sign_diffexp_RvP.tsv", col_select = 1) %>%
  pull(1)

PvE_GO_results <- enrichGO(gene = PvE_diffexp_genes,
                           OrgDb = org.Dmagna.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.01,
                           qvalueCutoff = 0.05,
                           keyType = "GID")
as.data.frame(PvE_GO_results)
PvE_fit <- plot(barplot(PvE_GO_results, showCategory = 20))
write.csv(PvE_GO_results, "PvE_GO_Results.csv")

RvE_GO_results <- enrichGO(gene = RvE_diffexp_genes,
                           OrgDb = org.Dmagna.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.01,
                           qvalueCutoff = 0.05,
                           keyType = "GID")
as.data.frame(RvE_GO_results)
RvE_fit <- plot(barplot(RvE_GO_results, showCategory = 20))
write.csv(RvE_GO_results, "RvE_GO_Results.csv")


#no results with 0.01 & 0.1 pvalue, so trying 0.2
RvP_GO_results <- enrichGO(gene = RvP_diffexp_genes,
                           OrgDb = org.Dmagna.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.2,
                           qvalueCutoff = 1,
                           keyType = "GID")
as.data.frame(RvP_GO_results)
RvP_fit <- plot(barplot(RvP_GO_results, showCategory = 20))
write.csv(RvP_GO_results, "RvP_GO_Results.csv")
