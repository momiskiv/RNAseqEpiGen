library(readr)
library(GenomicRanges)

setwd("~/final_project/github/Diff_Meth")

methylation_data <- read_csv("./Recovery_Pesticide/diff_meth_RvP.csv")

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

#Convert GFF data to GRanges object
gff_gr <- GRanges(
  seqnames = gff_data$seqid,
  ranges = IRanges(start = gff_data$start, end = gff_data$end),
  strand = gff_data$strand
)

# Find overlaps between methylation regions and GFF data
overlaps <- findOverlaps(methylation_gr, gff_gr)

# Extract overlapping features
overlapping_features <- gff_data[subjectHits(overlaps), ]

# Display overlapping features
print(overlapping_features)

write.csv(overlapping_features, file = "sign_diff_meth_features.csv")