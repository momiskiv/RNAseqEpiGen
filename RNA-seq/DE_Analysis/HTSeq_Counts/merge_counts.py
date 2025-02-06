## Short code to create a matrix of all the count files to optimise DESeqObject creation

import pandas as pd
import glob

# Loading count files from wd
count_files = glob.glob('*.counts')

sample_counts = {}

# Read each count file and store in the dict
for file in count_files:
    sample_name = file.replace('.counts', '')
    df = pd.read_csv(file, sep='\t', header=None, names=['Gene', sample_name])
    sample_counts[sample_name] = df.set_index('Gene')[sample_name]

# Combine series within df
if sample_counts:
    combined_df = pd.DataFrame(sample_counts)
    combined_df.insert(0, 'Gene', combined_df.index)
    combined_df.reset_index(drop=True, inplace=True)
    combined_df.to_csv('combined_counts_matrix.csv', index=False)
    print(combined_df.head())  # Display the first few rows of the combined dataframe
else:
    print('No .counts files found to combine.')
