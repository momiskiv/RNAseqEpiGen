library(readr)
library(GenomicRanges)
library(stringr)

setwd("~/final_project/github/Diff_Meth")

methylation_data <- read_csv("./Recovery_Eutrophic/diff_meth_RvE.csv")

gff_data <- read.table("NIES_genome.gff",
                       header = FALSE,
                       sep = "\t",
                       comment.char = "#")
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")


#=========================================================================================

# Calculate absolute methylation difference
methylation_data$abs_meth_diff <- abs(methylation_data$meth.diff)

# Sort data by absolute methylation difference in descending order
most_diff_meth_regions <- methylation_data[order(-methylation_data$abs_meth_diff), ]

# Display the top 10 regions
head(most_diff_meth_regions, 10)


#===============================================================================================

# Convert methylation data to GRanges object
methylation_gr <- GRanges(
  seqnames = methylation_data$chr,
  ranges = IRanges(start = methylation_data$start, end = methylation_data$end),
  strand = methylation_data$strand
)

#Filtering for just gene features
gff_genes <- gff_data[gff_data$type == "gene", ]

#Convert GFF data to GRanges object
gff_genes_gr <- GRanges(
  seqnames = gff_genes$seqid,
  ranges = IRanges(start = gff_genes$start, end = gff_genes$end),
  strand = gff_genes$strand
)

# Find overlaps between methylation regions and GFF data
overlaps <- findOverlaps(methylation_gr, gff_genes_gr)

# Extract overlapping features
overlapping_genes <- gff_genes[subjectHits(overlaps), ]

# Display overlapping genes
head(overlapping_genes)

write.csv(overlapping_genes, file = "RvE_sign_diff_meth_gene_features.csv")

#Extract names from attribute column
overlapping_genes$symbol <- str_extract(overlapping_genes$attributes, "(?<=Name=)[^;]+")

#Add to new df
diffmeth_genes <- data.frame(
  symbol = overlapping_genes$symbol,
  seqid = overlapping_genes$seqid)

write.csv(diffmeth_genes, file = "RvE_sign_diff_meth_genes.csv", row.names = FALSE)