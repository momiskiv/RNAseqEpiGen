#------------------------------------------------------------------
# Differential Gene Expression
#------------------------------------------------------------------

setwd("~/final_project/github/DE_Analysis")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Make sure to have the packages installed before loading!
#Some packages will require use of BiocManager 

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(apeglm)

#------------------------------------------------------------------
# Major support from tutorial:
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html

# Setting up data: need a matrix of counts and metadata
# Create matrix from HTSeq counts using python script on Git

cts <- as.matrix(read.csv("counts.csv", row.names="gene_id"))

coldata <- read.csv("Metadata.csv", row.names="Run")
coldata <- coldata[,c("Population","Biological_Replicate","Replicate")]
coldata$Population <- factor(coldata$Population)
coldata$Biological_Replicate <- factor(coldata$Biological_Replicate)
coldata$Replicate <- factor(coldata$Replicate)

# Order needs to be the same! IMPORTANT

rownames(coldata)
colnames(cts)
all(rownames(coldata) %in% colnames(cts))
#if FALSE, please ensure the sample names are the same
#if TRUE, continue
all (rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)] #reordered based on row order of metadata
all(rownames(coldata) == colnames(cts)) #should be TRUE; don't proceed if otherwise

# Construct of DESeq dataset

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Population)
#------------------------------------------------------------------

# Pre-filtering samples with low counts

nrow(dds) #26,140

dds = dds[rowMeans(counts(dds)) >10, ]
# 1,760/26,140
#USING THIS

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep, ]
#12,497/26,140; the DE was not so good with so many low-count genes

#Factor levels; R chooses reference levels based on alphabetical order

dds$Population <- factor(dds$Population, levels = c("Pesticide", "Recovery", "Eutrophic"))

# rlog transform counts (based on dds object design info)
rld = rlog(dds, blind=FALSE)

#------------------------------------------------------------------
# PCA plot

data = plotPCA(rld, intgroup = c("Population"), returnData=TRUE) #500 genes
percentVar = round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=Population)) + geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed(ratio=1,clip = "on")+
  geom_text_repel(aes(label=name), size=3,show.legend=FALSE, 
                  point.padding = 0.5, box.padding = 0.5,
                  segment.color = 'grey50', max.overlaps = Inf) +
  scale_colour_brewer(palette = "Set2") +
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.text = element_text(size=12))

# Scree plot

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))] # Top variable genes
pca <- prcomp(t(assay(rld)[select,]))
percentVar <- 100*(pca$sdev^2 / sum( pca$sdev^2 )) # the contribution to the total variance for each component

scree_plot <- data.frame(Variance = percentVar,
                         PC_Number = seq_along(percentVar))

ggplot(scree_plot, aes(x=PC_Number, y=Variance)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_bw() +
  ylab("% Variance") +
  xlab("Principal Component") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14))

# Top 3 PC = most relevant

# Heatmap of the top variable genes
my_colors = list(
  Population = c(Pesticide = "#66C2A5", Recovery = "#FC8D62", Eutrophic = "#8DA0CB"))

mat = assay(rld)[select, ]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(rld)[, "Population", drop=FALSE])
colnames(df2) <- c("Population")

pheatmap(mat,
         annotation_col=df2,
         show_rownames = FALSE,
         fontsize = 10,
         annotation_colors = my_colors,
         color = colorRampPalette(c("navy","white","firebrick3"))(50),
         main = "Top Variable Genes",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")


# two 1st samples plotted against each other to check consistency (for rlog and log2)

par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(1,2)] + 1),
     main = "log2 transformation", pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(1,2)],
     main = "rlog transformation", pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par( mfrow = c( 1, 1 ) )

#rlog better transformation

# Estimate size factors = normalize for dispersion

dds = DESeq2::estimateSizeFactors(dds)
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts",
             ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)

# Check sample distances - matrix only

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T,
         fontsize = 10)

# Check sample distances - Poisson distribution

poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
sample_names <- paste(rld$dex, rld$cell, sep = " - ")
colnames(samplePoisDistMatrix) <- colnames(sampleDistMatrix) #otherwise just shows numbers
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T,
         show_rownames = F,
         fontsize = 10)

#------------------------------------------------------------------

# Differential expression (used Benjamini-Hochberg adjustment;  p= <0.1)

dds = DESeq2::DESeq(dds, parallel=TRUE)
resultsNames(dds)

# RESULTS

# Recovery vs Pesticide
res_rec_v_pes <- results(dds, contrast = c("Population", "Recovery", "Pesticide"))
summary(res_rec_v_pes)
# out of 1760 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 206, 12%
# LFC < 0 (down)     : 142, 8.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 6)
#TOTAL DIFF.EXPR = 20.10%

# Recovery vs Eutrophic 
res_rec_v_eut <- results(dds, contrast = c("Population", "Recovery", "Eutrophic"))
summary(res_rec_v_eut)
# out of 1760 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 237, 13%
# LFC < 0 (down)     : 88, 5%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 6)
#TOTAL DIFF.EXPR = 18%

# Pesticide vs Eutrophic
res_pes_v_eut <- results(dds, contrast = c("Population", "Pesticide", "Eutrophic"))
summary(res_pes_v_eut)
# out of 1760 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 373, 21%
# LFC < 0 (down)     : 301, 17%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 6)
#TOTAL DIFF.EXPR = 38%

# Distribution of coefficents of the model
# MA plots

# Recovery vs Pesticide
plotMA(res_rec_v_pes, alpha=0.1, ylim=c(-8,8),cex=1.2, cex.lab=1.5, cex.axis=1.5,
       xlab='Mean of Normalised Counts', ylab='Log2 Fold Change', main = "Recovery vs Pesticide (p=<0.1)",
       colSig="blue",colNonSign="gray50")
#makeup
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
abline(h=0, col="black", lwd=2, lty=2)
legend("topright", legend = c("Upregulated","Non-significant"),
       col =c("blue","gray50"),
       pch =16, pt.cex = 1, cex = 0.8, bg = "white")

# Recovery vs Eutrophic 
plotMA(res_rec_v_eut, alpha=0.1, ylim=c(-5,5),cex=1.2, cex.lab=1.5, cex.axis=1.5,
       xlab='Mean of Normalised Counts', ylab='Log2 Fold Change', main = "Recovery vs Eutrophic (p=<0.1)",
       colSig="blue",colNonSign="gray50")
#makeup
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
abline(h=0, col="black", lwd=2, lty=2)
legend("topright", legend = c("Upregulated","Non-significant"),
       col =c("blue","gray50"),
       pch =16, pt.cex = 1, cex = 0.8, bg = "white")

# Pesticide vs Eutrophic
plotMA(res_pes_v_eut, alpha=0.1, ylim=c(-10,10),cex=1.2, cex.lab=1.5, cex.axis=1.5,
       xlab='Mean of Normalised Counts', ylab='Log2 Fold Change', main = "Pesticide vs Eutrophic (p=<0.1)",
       colSig="blue",colNonSign="gray50")
#makeup
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
abline(h=0, col="black", lwd=2, lty=2)
legend("topright", legend = c("Upregulated","Non-significant"),
       col =c("blue","gray50"),
       pch =16, pt.cex = 1, cex = 0.8, bg = "white")

# plot of p-vals excluding genes with very small counts
#HISTOGRAMS

# Recovery vs Pesticide 
hist(res_rec_v_pes$pvalue[res_rec_v_pes$baseMean > 1 & res_rec_v_pes$pvalue <= 1], 
     breaks = seq(0, 1, by = 0.05),
     col = "darkblue",
     border = "white",
     xlab = "P-value", 
     ylab = "Frequency",
     cex.axis = 1.45, 
     cex.lab = 1.45,
     main = "Recovery vs Pesticide")

# Recovery vs Eutrophic 
hist(res_rec_v_eut$pvalue[res_rec_v_eut$baseMean > 1 & res_rec_v_eut$pvalue <= 1], 
     breaks = seq(0, 1, by = 0.05),
     col = "darkblue",
     border = "white",
     xlab = "P-value", 
     ylab = "Frequency",
     cex.axis = 1.45, 
     cex.lab = 1.45,
     main = "Recovery vs Eutrophic")

# Pesticide vs Eutrophic
hist(res_pes_v_eut$pvalue[res_pes_v_eut$baseMean > 1 & res_pes_v_eut$pvalue <= 1], 
     breaks = seq(0, 1, by = 0.05),
     col = "darkblue",
     border = "white",
     xlab = "P-value", 
     ylab = "Frequency",
     cex.axis = 1.45, 
     cex.lab = 1.45,
     main = "Pesticide vs Eutrophic")

# Order the results by fold change to see the biggest changing genes

#RvP
res_RvP_ordered=res_rec_v_pes[order(res_rec_v_pes$log2FoldChange),]
head(res_RvP_ordered)
res_RvP_ordered<-as.data.frame(res_RvP_ordered)
res_RvP_ordered$gene<-rownames(res_RvP_ordered)
rownames(res_RvP_ordered)<-c()
write.table(as.data.frame(res_RvP_ordered), file="diff_exp_output_log2FC_all_genes_RvP.tsv",
            sep="\t", quote = F, col.names = T, row.names = F)

#PvE
res_PvE_ordered=res_pes_v_eut[order(res_pes_v_eut$log2FoldChange),]
head(res_PvE_ordered)
res_PvE_ordered<-as.data.frame(res_PvE_ordered)
res_PvE_ordered$gene<-rownames(res_PvE_ordered)
rownames(res_PvE_ordered)<-c()
write.table(as.data.frame(res_PvE_ordered), file="diff_exp_output_log2FC_all_genes_PvE.tsv",
            sep="\t", quote = F, col.names = T, row.names = F)

#RvE
res_RvE_ordered=res_rec_v_eut[order(res_rec_v_eut$log2FoldChange),]
head(res_RvE_ordered)
res_RvE_ordered<-as.data.frame(res_RvE_ordered)
res_RvE_ordered$gene<-rownames(res_RvE_ordered)
rownames(res_RvE_ordered)<-c()
write.table(as.data.frame(res_RvE_ordered),
            file="diff_exp_output_log2FC_all_genes_RvE.tsv",
            sep="\t", quote = F, col.names = T, row.names = F)

# out of 1760 with nonzero total read count

#RvE
res_RvE_significant<- subset(res_rec_v_eut, log2FoldChange > 1.5 | log2FoldChange < -1.5)
res_RvE_significant <- subset(res_RvE_significant, padj < 0.05)
nrow(res_RvE_significant) # 68 total differentially expressed genes
nrow(res_RvE_significant[res_RvE_significant$log2FoldChange > 1.5,]) #64 upregulated in recovery
nrow(res_RvE_significant[res_RvE_significant$log2FoldChange < -1.5,]) #4 upregulated in eutrophic

write.table(as.data.frame(res_RvE_significant),
            file="sign_diffexp_RvE.tsv",
            sep="\t", quote = F, col.names = NA, row.names=T)

#RvP
res_RvP_significant<- subset(res_rec_v_pes, log2FoldChange > 1.5 | log2FoldChange < -1.5)
res_RvP_significant <- subset(res_RvP_significant, padj < 0.05)
nrow(res_RvP_significant) # 78 total differentially expressed genes
nrow(res_RvP_significant[res_RvP_significant$log2FoldChange > 1.5,]) #53 upregulated in recovery
nrow(res_RvP_significant[res_RvP_significant$log2FoldChange < -1.5,]) #25 upregulated in eutrophic

write.table(as.data.frame(res_RvE_significant),
            file="sign_diffexp_RvP.tsv",
            sep="\t", quote = F, col.names = NA, row.names=T)

#PvE
res_PvE_significant<- subset(res_pes_v_eut, log2FoldChange > 1.5 | log2FoldChange < -1.5)
res_PvE_significant <- subset(res_PvE_significant, padj < 0.05)
nrow(res_PvE_significant) # 150 total differentially expressed genes
nrow(res_PvE_significant[res_PvE_significant$log2FoldChange > 1.5,]) # 90 upregulated in recovery
nrow(res_PvE_significant[res_PvE_significant$log2FoldChange < -1.5,]) # 60 upregulated in eutrophic

write.table(as.data.frame(res_RvE_significant),
            file="sign_diffexp_PvE.tsv",
            sep="\t", quote = F, col.names = NA, row.names=T)

#------------------------------------------------------------------

# Heatmap of top differentially expressed

#RvE
my_colors = list(
  Population = c(Recovery = "#FC8D62", Eutrophic = "#8DA0CB"))

df2 = as.data.frame(colData(rld)[,c("Population"),drop=FALSE])
df2 = df2[df2$Population %in% c("Recovery","Eutrophic"), ,drop=FALSE]

n=50
topdiff = head(c(1:nrow(res_rec_v_eut))[order(res_rec_v_eut$padj)],n)
mat = assay(rld)[ topdiff, ]
mat = mat[,colnames(mat) %in% rownames(df2)]
mat = mat - rowMeans(mat)

pheatmap(mat, 
         annotation_col = df2,
         show_rownames = F,
         fontsize = 12,
         annotation_colors = my_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Recovery vs Eutrophic (Top 50)",
         annotation_names_col = F,
         annotation_legend = T,
         legend = TRUE)

#PvE
my_colors = list(
  Population = c(Pesticide = "#66C2A5",Eutrophic = "#8DA0CB"))

df2 = as.data.frame(colData(rld)[,c("Population"),drop=FALSE])
df2 = df2[df2$Population %in% c("Pesticide","Eutrophic"), ,drop=FALSE]

n=50
topdiff = head(c(1:nrow(res_pes_v_eut))[order(res_pes_v_eut$padj)],n)
mat = assay(rld)[ topdiff, ]
mat = mat[,colnames(mat) %in% rownames(df2)]
mat = mat - rowMeans(mat)

pheatmap(mat, 
         annotation_col = df2,
         show_rownames = F,
         fontsize = 12,
         annotation_colors = my_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Pesticide vs Eutrophic (Top 50)",
         annotation_names_col = F,
         annotation_legend = T,
         legend = TRUE)

#RvP
my_colors = list(
  Population = c(Pesticide = "#66C2A5", Recovery = "#FC8D62"))

df2 = as.data.frame(colData(rld)[,c("Population"),drop=FALSE])
df2 = df2[df2$Population %in% c("Recovery","Pesticide"), ,drop=FALSE]

n=50
topdiff = head(c(1:nrow(res_rec_v_pes))[order(res_rec_v_pes$padj)],n)
mat = assay(rld)[ topdiff, ]
mat = mat[,colnames(mat) %in% rownames(df2)]
mat = mat - rowMeans(mat)

pheatmap(mat, 
         annotation_col = df2,
         show_rownames = F,
         fontsize = 12,
         annotation_colors = my_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Recovery vs Pesticide (Top 50)",
         annotation_names_col = F,
         annotation_legend = T,
         legend = TRUE)

#------------------------------------------------------------------

# Plot diff expressed

#PvE
n=15 # 150 diff.expr. genes -> 10% = 15
selGenes = head(rownames(res_pes_v_eut)[order(res_pes_v_eut$padj)],n)
data = do.call(rbind,lapply(selGenes, function(gene) data.frame(gene=gene,plotCounts(dds, gene=gene, intgroup=c("Population"), returnData=TRUE))))
data = data[data$Population %in% c("Pesticide","Eutrophic"), ]

ggplot(data, aes(x=Population, y=count, fill=Population)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Population") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=14, face="bold"),
        legend.text=element_text(size=12),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=16, hjust=0.5),)+
  scale_fill_manual("Population", 
                    breaks=c("Pesticide", "Eutrophic"),
                    values=c("#66C2A5", "#8DA0CB")) +
  ggtitle("Pesticide vs Eutrophic (Top 10%)")


#RvE
n=7 # 68 diff.expr. genes -> 10% = 6.8 -> 7 (1 s.f.)
selGenes = head(rownames(res_rec_v_eut)[order(res_rec_v_eut$padj)],n)
data = do.call(rbind,lapply(selGenes, function(gene) data.frame(gene=gene,plotCounts(dds, gene=gene, intgroup=c("Population"), returnData=TRUE))))
data = data[data$Population %in% c("Recovery","Eutrophic"), ]

ggplot(data, aes(x=Population, y=count, fill=Population)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Population") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=14, face="bold"),
        legend.text=element_text(size=12),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=16, hjust=0.5),)+
  scale_fill_manual("Population", 
                    breaks=c("Recovery", "Eutrophic"),
                    values=c("#FC8D62", "#8DA0CB")) +
  ggtitle("Recovery vs Eutrophic (Top 10%)")

#RvP
n=8 # 78 diff.expr. genes -> 10% = 7.8 -> 8 (1 s.f.)
selGenes = head(rownames(res_rec_v_pes)[order(res_rec_v_pes$padj)],n)
data = do.call(rbind,lapply(selGenes, function(gene) data.frame(gene=gene,plotCounts(dds, gene=gene, intgroup=c("Population"), returnData=TRUE))))
data = data[data$Population %in% c("Recovery","Pesticide"), ]

ggplot(data, aes(x=Population, y=count, fill=Population)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Population") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=14, face="bold"),
        legend.text=element_text(size=12),
        axis.text.y=element_text(size=10),
        plot.title = element_text(size=16, hjust=0.5),)+
  scale_fill_manual("Population", 
                    breaks=c("Recovery", "Pesticide"),
                    values=c("#FC8D62", "#66C2A5")) +
  ggtitle("Recovery vs Pesticide (Top 10%)")

#------------------------------------------------------------------

# Volcano plot

#RvE
res_RvE_df <- as.data.frame(res_rec_v_eut)
res_RvE_df $gene <- row.names(res_RvE_df)

res_RvE_df$expression <- "Insignificant"
res_RvE_df$expression[res_RvE_df$log2FoldChange > 1.5 & res_RvE_df$padj < 0.05] <- "Upregulated"
res_RvE_df$expression[res_RvE_df$log2FoldChange < -1.5 & res_RvE_df$padj < 0.05] <- "Downregulated"

# Top 7 genes marked
top_genes <- res_RvE_df[order(res_RvE_df$padj), ][1:7, ]

ggplot(res_RvE_df, aes(x = log2FoldChange, y = -log10(padj), color = expression, show.legend = T)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(
    values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey"), name = "Gene expression") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "grey40") +
  labs(
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"(adjusted p-value)"),
    title = "Recovery vs Eutrophic") +
  guides(color = guide_legend(override.aes = list(label = "")))
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))

#RvP
  res_RvP_df <- as.data.frame(res_rec_v_pes)
  res_RvP_df $gene <- row.names(res_RvP_df)
  
  res_RvP_df$expression <- "Insignificant"
  res_RvP_df$expression[res_RvP_df$log2FoldChange > 1.5 & res_RvP_df$padj < 0.05] <- "Upregulated"
  res_RvP_df$expression[res_RvP_df$log2FoldChange < -1.5 & res_RvP_df$padj < 0.05] <- "Downregulated"
  
  # Top 8 genes marked
  top_genes <- res_RvP_df[order(res_RvP_df$padj), ][1:8, ]
  
  ggplot(res_RvP_df, aes(x = log2FoldChange, y = -log10(padj), color = expression, show.legend = T)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey"), name = "Gene expression") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      max.overlaps = Inf,
      box.padding = 0.5,
      segment.color = "grey40") +
    labs(
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(adjusted p-value)"),
      title = "Recovery vs Pesticide") +
    guides(color = guide_legend(override.aes = list(label = "")))
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12))

#PvE
  res_PvE_df <- as.data.frame(res_pes_v_eut)
  res_PvE_df$gene <- row.names(res_PvE_df)
  
  res_PvE_df$expression <- "Insignificant"
  res_PvE_df$expression[res_PvE_df$log2FoldChange > 1.5 & res_PvE_df$padj < 0.05] <- "Upregulated"
  res_PvE_df$expression[res_PvE_df$log2FoldChange < -1.5 & res_PvE_df$padj < 0.05] <- "Downregulated"
  
  # Top 15 genes marked
  top_genes <- res_PvE_df[order(res_PvE_df$padj), ][1:15, ]
  
  ggplot(res_PvE_df, aes(x = log2FoldChange, y = -log10(padj), color = expression, show.legend = T)) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey"), name = "Gene expression") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      max.overlaps = Inf,
      box.padding = 0.5,
      segment.color = "grey40") +
    labs(
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(adjusted p-value)"),
      title = "Pesticide vs Eutrophic") +
    guides(color = guide_legend(override.aes = list(label = "")))
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12))

#------------------------------------------------------------------
# Pull out a list of just the differentially expressed genes for the supplementary
# Use Python script -> saves as .txt (tab separated)


