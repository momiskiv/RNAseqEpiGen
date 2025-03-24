
######################################
###### ML FEEDING DATA CREATION ######
######################################

import pandas as pd

# Load data
perc_meth = pd.read_csv("~/final_project/github/Diff_Meth/perc_meth.csv", index_col=0)
geno_env = pd.read_excel("~/final_project/Data/genotypes_and_environmental_data.xlsx", sheet_name="genotype_dates")

# Map sample names
sample_mapping = {
    'recovery_0_1':'CWP_LRV0_1',
    'recovery_0_2':'CWP_LRV0_2',
    'recovery_0_4':'CWP_LRV0_4',
    'recovery_1_2':'CWP_LRV1_2',
    'recovery_2_1':'CWP_LRV2_1',
    'recovery_2_5_9':'CWP_LRV2_5_9',
    'recovery_2_5_11':'CWP_LRV2_5_11',
    'recovery_3_5_1':'CWP_LRV3_5_1',
    'recovery_3_5_2':'CWP_LRV3_5_2',
    'recovery_3_5_15':'CWP_LRV3_5_15',
    'recovery_3_6':'CWP_LRV3_6',
    'pesticide_6_2':'PP_LRV6_2',
    'pesticide_6_3':'PP_LRV6_3',
    'pesticide_7_3':'PP_LRV7_3',
    'pesticide_7_5':'PP_LRV7_5',
    'pesticide_7_5_4':'PP_LRV7_5_4',
    'pesticide_8_5_3':'PP_LRV8_5_3',
    'pesticide_9_5_1':'PP_LRV9_5_1',
    'pesticide_9_5_3':'PP_LRV9_5_3',
    'pesticide_9_6':'PP_LRV9_6',
    'pesticide_9_20':'PP_LRV9_20',
    'eutrophic_12_3':'EP_LRV12_3',
    'eutrophic_12_4':'EP_LRV12_4',
    'eutrophic_12_5_1':'EP_LRV12_5_1',
    'eutrophic_13_1':'EP_LRV13_1',
    'eutrophic_13_2':'EP_LRV13_2',
    'eutrophic_13_3':'EP_LRV13_3',
    'eutrophic_13_5_1':'EP_LRV13_5_1',
    'eutrophic_14_5_1':'EP_LRV14_5_1',
    'eutrophic_15_5_1':'EP_LRV15_5_1'
}

# Choose pollutant (fungicide, herbicide, insecticide, pesticide)
pollutant = "pesticide_lvl"

# Convert geno_env into a dict
pollutant_levels = geno_env.set_index("full_sample_name")[pollutant].to_dict()
population_data = geno_env.set_index("full_sample_name")['population'].to_dict()

# Select only relevant columns from perc_meth
perc_meth_filtered = perc_meth.loc[:, sample_mapping.keys()].T 

# Map sample names and get pollutant levels
perc_meth_filtered.index = perc_meth_filtered.index.map(sample_mapping)

# Add pollutant levels
perc_meth_filtered["pesticide"] = perc_meth_filtered.index.map(pollutant_levels)
perc_meth_filtered["population"] = perc_meth_filtered.index.map(population_data)

# Reorder columns to place pollutant first
final_df = perc_meth_filtered.reset_index(drop=True)
cols = ["pesticide"] + [col for col in final_df.columns if col not in ["pesticide", "population"]] + ["population"]
final_df = final_df[cols]

# Save to CSV
final_df.to_csv(f"{pollutant}_ml_training_data.csv", index=False)