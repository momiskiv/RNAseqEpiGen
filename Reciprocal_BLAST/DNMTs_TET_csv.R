## Getting counts for reciprocal BLAST searches for DNMTs and TET

setwd("/home/nm471/final_project/github/Reciprocal_BLAST/")
library(readr)

# Genes
genes <- c("dnmt1", "dnmt3a", "tet")

for (gene in genes) {

  file_name <- paste0("dmagna_", gene, "_reciprocal.txt")

  data <- read_delim(file_name, "\t", escape_double = F, trim_ws = T, col_names = F)
  
  results <- table(data$X2)
  results_df <- as.data.frame(results)
  colnames(results_df) <- c("Refseq", "Count")
  results_df <- results_df[order(-results_df$Count), ] #descending order
  
  cat(gene, ":\n")
  print(head(results_df, 3)) #top 3 hits
  
  #Annotation columns
  results_df$Gene_name <- ""
  results_df$Gene_symbol <- ""
  results_df$Gene_description <- ""
  
  #CSV
  output_f <- paste0("dmagna_",gene,"_matches.csv")
  write.csv(head(results_df, 3), output_f, row.names = F)
}