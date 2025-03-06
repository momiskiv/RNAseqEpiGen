
## -------------------------------------------------------------------------
## All Populations - ML 
## -------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE, 
               repos=BiocManager::repositories(),
               dependencies=TRUE)

library(grid) # rewrite of the graphics layout capabilities, plus some support for interaction.
library(readr) 
library(methylKit) # finding differentially methylated genes, similar to Deseq2
library(ggfortify) # plotting tools for statistics
library(future.apply)
library(future)

setwd("~/final_project/github/Diff_Meth")

## -------------------------------------------------------------------------
# GET METHYLOBJ FOR SAMPLES# 

# All in one folder; changed from specific diff_meth
sample.list <- as.list(list.files(
  path = "./Coverage_Files",
  pattern = ".*\\.CpG_report\\.merged_CpG_evidence\\.cov$",
  full.names = T))


# In treatment, specify how many representatives in each groups compared (group 0 and 1)
# No filtering here, done in following step
# Includes all populations
CPGRaw <- methRead(sample.list, 
                   sample.id = list("recovery_0_1",
                                    "recovery_0_2",
                                    "recovery_0_4",
                                    "recovery_1_2",
                                    "recovery_2_1",
                                    "recovery_2_5_9",
                                    "recovery_2_5_11",
                                    "recovery_3_5_1",
                                    "recovery_3_5_2",
                                    "recovery_3_5_15",
                                    "recovery_3_6",
                                    "pesticide_6_2",
                                    "pesticide_6_3",
                                    "pesticide_7_3",
                                    "pesticide_7_5",
                                    "pesticide_7_5_4",
                                    "pesticide_8_5_3",
                                    "pesticide_9_5_1",
                                    "pesticide_9_5_3",
                                    "pesticide_9_6",
                                    "pesticide_9_20",
                                    "eutrophic_12_3",
                                    "eutrophic_12_4",
                                    "eutrophic_12_5_1",
                                    "eutrophic_13_1",
                                    "eutrophic_13_2",
                                    "eutrophic_13_3",
                                    "eutrophic_13_5_1",
                                    "eutrophic_14_5_1",
                                    "eutrophic_15_5_1"
                   ),
                   assembly="D.magna",
                   treatment=c(rep(0,11), rep(1,10), rep(2,9)),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= getwd())

gc()

# DATA FORMATTING # 

# Filter by coverage
filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

# Select only CpGs found in all samples
meth_all_data <- unite(filtered_data, destrand=F) 
nrow(meth_all_data)
# 120,863

gc()

# Pull data into a usable dataframe
df_meth_all <- getData(meth_all_data)

# remove first few uninformative columns
df_meth_all <- df_meth_all[,-c(1:4)]

# subset data into individual dataframes per sample
all_data <- lapply(seq(1, ncol(df_meth_all), by=3), function(i) 
  df_meth_all[i: pmin((i+1), ncol(df_meth_all))])

# STATISTICAL TESTING #

# Define binomial test
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}

# Set up output dataframe
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

# Parallelisation plan
plan(multisession, workers = 4) #change depending on system specs; availableCores()

gc() # make some space before parallel loop

# Running binomial on every row parallely
binomial_data <- future_lapply(seq_along(all_data), function(i) {
  data <- all_data[[i]]
  data <- data[complete.cases(data),] # only rows with values
  colnames(data) <- c("CT", "Ccount")
  data$pVal <- mapply(bt, data$Ccount, data$CT)
  data$FDR <- p.adjust(data$pVal, method = "BH", n = length(data$pVal)) # false discovery rate adj
  data$row <- as.numeric(rownames(data))
  dfmeth <- subset(data, data$FDR < 0.05)
  return(dfmeth)
})
# 00:01:21:40

gc()

# Combine the results from all iterations
allrows <- do.call(rbind, binomial_data)

# Get positions which are methylated in at least one sample
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) # 3,002

# Keep only these positions for diff meth test
subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

gc()

## S4 method for signature 'methylBase'
perc_meth <- percMethylation(subset_methBase, rowids = TRUE)
head(perc_meth)
write.csv(perc_meth, file = "perc_meth.csv")