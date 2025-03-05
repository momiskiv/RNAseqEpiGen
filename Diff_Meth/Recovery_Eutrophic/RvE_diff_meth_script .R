
## -------------------------------------------------------------------------
## Recovery v Eutrophic
## -------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(devtools)
install_github("al2na/methylKit", build_vignettes=FALSE, 
               repos=BiocManager::repositories(),
               dependencies=TRUE)

library(grid) # rewrite of the graphics layout capabilities, plus some support for interaction.
library(readr) 
library(ggplot2) 
library(methylKit) # finding differentially methylated genes, similar to Deseq2
library(ggdendro) # set of tools for dendrograms and tree plots using 'ggplot2'
library(ggfortify) # plotting tools for statistics
library(ggrepel)
library(viridis) # colour palettes
library(future.apply)
library(future)

setwd("~/final_project/github/Diff_Meth")

## -------------------------------------------------------------------------
# GET METHYLOBJ FOR SAMPLES# 

# Using regex; change pattern depending on own file names
sample.list <- as.list(list.files(
  path = "./Coverage_Files/RvE",
  pattern = ".*\\.CpG_report\\.merged_CpG_evidence\\.cov$",
  full.names = T))


# In treatment, specify how many representatives in each groups compared (group 0 and 1)
# No filtering here, done in following step
# Edit with own sample names
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
                   treatment=c(rep(0,11), rep(1,9)),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir= getwd())

# DATA FORMATTING # 

# Filter by coverage
filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)

# Select only CpGs found in all samples
meth_all_data <- unite(filtered_data, destrand=F) 
nrow(meth_all_data) 
# 148,395

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

# Combine the results from all iterations
allrows <- do.call(rbind, binomial_data)

# Get positions which are methylated in at least one sample
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) # 4,302

# Keep only these positions for diff meth test
subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

# Write out data for later use
subset_methBase_out <- getData(subset_methBase)
write.table(subset_methBase_out, file="RvE_sites_for_PCA.txt", sep="\t", quote=F, row.names = F, col.names = T)

## -------------------------------------------------------------------------
# PCA PLOT #

PCA_data <- PCASamples(subset_methBase, screeplot=F, obj.return = T)
PCA_data1 <- as.data.frame(PCA_data$x)
PCA_data1$sample <- row.names(PCA_data1)

#Edit name and representatives based on comparison
PCA_data1$Condition <- c(rep("Recovery", 11), rep("Eutrophic", 9))

percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )

#Edit title & Legend
ggplot(PCA_data1, aes(PC1, PC2, colour=Condition))+
  geom_point(size=10, alpha = 0.8) +
  scale_color_viridis(discrete = T) +
  theme_bw()+
  xlab(paste0("PC1:",percentage[1],"variance")) +
  ylab(paste0("PC2:",percentage[2],"variance")) +
  theme(axis.text=element_text(size=26),
        axis.title=element_text(size=28),
        legend.text=element_text(size=30),
        legend.title=element_blank(),
        title=element_text(size=26))+
  ggtitle("Recovery v Eutrophic")

# SCREE PLOT #

eigenvalues <- PCA_data$sdev^2
percentage <- round(eigenvalues / sum(eigenvalues) * 100, 1)

df <- data.frame(PC = factor(seq_along(eigenvalues)), Eigenvalue = eigenvalues)

ggplot(df, aes(x = PC, y = Eigenvalue)) +
  geom_bar(stat = "identity", fill = "skyblue", alpha = 0.7) +  
  geom_line(aes(group = 1), color = "blue", size = 1) +        
  geom_point(color = "blue", size = 2) +                   
  labs(title = "RvE Scree Plot",
       x = "Component Number", y = "Eigenvalue") +
  theme_minimal()

## -------------------------------------------------------------------------
# DIFFERENTIAL METHYLATION OF CpG #

diff_meth <- calculateDiffMeth(subset_methBase, method='qvalue')
diff_meth_filtered <- getMethylDiff(diff_meth, difference=15, qvalue=0.01) # can mess around with these values
nrow(diff_meth_filtered)
# 173

head(diff_meth_filtered)

write.csv(diff_meth_filtered, file="diff_meth_RvE.csv") #edit for comparison; keep format