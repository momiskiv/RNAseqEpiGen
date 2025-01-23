## Getting counts for reciprocal BLAST searches for DNMTs and TET

library(readr)

# Genes
genes <- c("DNMT1", "DNMT3", "TET")

for (gene in genes) {

  file_name <- paste0("dmagna_", gene, "_reciprocal.txt")

  data <- read_delim(file_name, "\t", escape_double = F, trim_ws = T, col_names = F)
  
  results <- table(data$X2)
  cat(gene, ":\n")
  print(results)
}