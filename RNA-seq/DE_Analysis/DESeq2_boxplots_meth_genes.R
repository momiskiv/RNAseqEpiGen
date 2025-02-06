#------------------------------------------------------------------
# Differential Gene Expression for Methylation-related genes of interest
#------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Make sure to have the packages installed before loading!
#Some packages will require use of BiocManager 

library(DESeq2)
library(tidyverse)
library(ggplot2)

#------------------------------------------------------------------
# Major support from tutorial:
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html

#Script is also identical to original DESeq2 R script, but I used different variable names as it's within the same project

# Setting up data: need a matrix of counts and metadata
# Create matrix from HTSeq counts using python script on Git

countData <- as.matrix(read.csv("counts.csv", row.names="gene_id"))

colData <- read.csv("Metadata.csv", row.names="Run")
colData <- colData[,c("Population","Biological_Replicate","Replicate")]
colData$Population <- factor(coldata$Population)
colData$Biological_Replicate <- factor(colData$Biological_Replicate)
colData$Replicate <- factor(colData$Replicate)

# Order needs to be the same! IMPORTANT

rownames(colData)
colnames(countData)
all(rownames(colData) %in% colnames(countData))
#if FALSE, please ensure the sample names are the same
#if TRUE, continue
all (rownames(colData) == colnames(countData))
countData <- countData[, rownames(colData)] #reordered based on row order of metadata
all(rownames(colData) == colnames(countData)) #should be TRUE; don't proceed if otherwise

# Construct of DESeq dataset

dds_v2 <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Population)

dds_v2 <- DESeq(dds_v2)
normalised_counts <- counts(dds_v2, normalized=TRUE)

#NO FILTERING! The counts for the methylation genes are very low, hence why there was no sight of them in the initial run
#Here we just want to focus on seeing how the expression, even if little, compares between the temporal populations

results_RP_v2 <- results(dds_v2, contrast = c("Population", "Recovery", "Pesticide"))
summary(results_RP_v2)

# out of 17134 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 348, 2%
# LFC < 0 (down)     : 169, 0.99%
# outliers [1]       : 0, 0%
# low counts [2]     : 10013, 58%
# (mean count < 1)

results_RE_v2 <- results(dds_v2, contrast = c("Population", "Recovery", "Eutrophic"))
summary(results_RE_v2)

# out of 17134 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 389, 2.3%
# LFC < 0 (down)     : 89, 0.52%
# outliers [1]       : 0, 0%
# low counts [2]     : 11628, 68%
# (mean count < 1)

results_PE_v2 <- results(dds_v2, contrast = c("Population", "Pesticide", "Eutrophic"))
summary(results_PE_v2)

# out of 17134 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 649, 3.8%
# LFC < 0 (down)     : 459, 2.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 10336, 60%
# (mean count < 1)

genes_of_interest <- c("LOC116925494", "LOC116922472", "LOC116923817", "LOC116917382")

# LOC116925494 = DNMT1
# LOC116922472 and LOC116923817 = DNMT3B
# LOC116917382 =  TET


for (gene in genes_of_interest) {
  gene_counts <- normalised_counts[gene, ]
  df <- data.frame(sample = names(gene_counts),
                   counts = gene_counts,
                   Population = colData$Population)
  p <- ggplot(df, aes(x = Population, y = counts, fill = Population)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold"), legend.position = "none") +
    labs(title = paste0("Expression of ", gene),
         y = "Normalized Counts",
         x = "Population") 
  
    ggsave(filename = paste0("gen_expr_",gene,".png"), plot = p, width = 10, height = 7, dpi = 300)
}