### simple parser for generated tsv from DESeq2 RNA-seq analysis to get list of differentially expressed genes with p-values <= 0.05

import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="DESeq2 tsv file processer to get a txt file of differentially expressed genes from RNA-seq with significant p-values.")
parser.add_argument("input_file", help = "Path to the input file. Has to be a tsv file from DESeq2 (columns gene, log2FoldChange, padj).")
parser.add_argument("-o", "--output_file", help = "Path to the output file. You can create it here, so choose any name you'd like!")
parser.add_argument("-v", "--verbose", action ="store_true", help = "Enable verbose mode. Allows to echo what will be written in the files within the command-line.")

args = parser.parse_args()

def process_deseq_tsv(input_file, output_file, verbose=False):

    df = pd.read_csv(input_file,sep="\t")
    sign_genes = df[df["padj"] <= 0.05][["gene", "log2FoldChange","padj"]]

    with open(output_file, "w") as f:
        for _,row in sign_genes.iterrows():
            gene = f"{row["gene"]}\t{row["log2FoldChange"]}\tp = {row["padj"]}\n"
            f.write(gene)
            if verbose:
                print(gene.strip())

input_file = args.input_file
output_file = args.output_file
verbose = args.verbose

process_deseq_tsv(input_file, output_file, verbose)
if verbose:
    print(f"Parsing complete. Output written to {output_file} in working directory.")