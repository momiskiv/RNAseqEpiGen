##short code initially made to check if the genes of interests were in the original count file from HTSeq

import pandas as pd

genes_of_interest = ["LOC116925494", "LOC116922472", "LOC116923817", "LOC116917382"]


df = pd.read_csv("counts.csv", index_col=0)

for gene in genes_of_interest:
    if gene in df.index:
        print(f"Gene {gene} found:")
        print(df.loc[gene])
        print("\n")
    else:
        print(f"Gene {gene} not found in the dataset.")
        print("\n")
