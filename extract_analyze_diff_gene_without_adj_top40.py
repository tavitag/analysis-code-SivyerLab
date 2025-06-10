# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 16:21:46 2023

@author: garretav
"""

import os
print(os.getcwd())
import pandas as pd
import csv

# Load RGC Atlas coordinates into a pandas DataFrame and set the first column as the index
RGC_Atlas_coordinates = pd.read_csv("RGC_Atlas_coordinates.txt", delimiter="\t", index_col=0)

# View the DataFrame (display the first few rows in Jupyter notebooks)
print(RGC_Atlas_coordinates.head())

# Load RGC Atlas into a pandas DataFrame (may take several  minutes)
RGC_Atlas = pd.read_csv("RGC_Atlas.csv")

# View the DataFrame 
print(RGC_Atlas.head())
#%% 
# Extract unique IDs for cells in a specific cluster
Clust_A = "44_Novel"
Clust_B = "18_Novel"

# Find all cells in the coordinate atlas part of the target clusters
Clust_A_cells = RGC_Atlas_coordinates[RGC_Atlas_coordinates["Cluster"] == Clust_A].index
Clust_B_cells = RGC_Atlas_coordinates[RGC_Atlas_coordinates["Cluster"] == Clust_B].index

Clust_A_uni_IDs = RGC_Atlas.columns[RGC_Atlas.columns.isin(Clust_A_cells)]
Clust_B_uni_IDs = RGC_Atlas.columns[RGC_Atlas.columns.isin(Clust_B_cells)]

Clust_A_filt = RGC_Atlas[["GENE"] + list(Clust_A_cells)]
Clust_B_filt = RGC_Atlas[["GENE"] + list(Clust_B_cells)]

# Save your filtered data to CSV files with cluster names in the file names
Clust_A_file_name = f"{Clust_A}_filt.csv"
Clust_B_file_name = f"{Clust_B}_filt.csv"

Clust_A_filt.to_csv(Clust_A_file_name, index=False, quoting=csv.QUOTE_NONNUMERIC)
Clust_B_filt.to_csv(Clust_B_file_name, index=False, quoting=csv.QUOTE_NONNUMERIC)

# Load the datasets 
df1 = pd.read_csv(f"{Clust_A_file_name}")
df2 = pd.read_csv(f"{Clust_B_file_name}")

# Display the first few rows of each dataset to confirm data structure
df1.head(), df2.head()

# Get the number of genes and samples in each dataset
num_genes_df1 = df1.shape[0]
num_samples_df1 = df1.shape[1] - 1  # subtract 1 for the 'GENE' column
num_genes_df2 = df2.shape[0]
num_samples_df2 = df2.shape[1] - 1  # subtract 1 for the 'GENE' column

num_genes_df1, num_samples_df1, num_genes_df2, num_samples_df2

# Check if the genes are the same and in the same order in both datasets. Should be okay if you used 
all(df1["GENE"] == df2["GENE"])
all(df1.iloc[:, 0] == df2.iloc[:, 0])

#load required packages
from scipy.stats import ttest_ind
import numpy as np

# Remove the 'GENE' column to only keep the expression data
data1 = df1.iloc[:, 1:]
data2 = df2.iloc[:, 1:]

# Transpose the data so that genes are columns and samples are rows
data1 = data1.transpose()
data2 = data2.transpose()

# Perform a t-test for each gene
t_statistics, p_values = ttest_ind(data1, data2, axis=0, equal_var=False, nan_policy='omit')

# Create a DataFrame with the t-statistics and p-values
diff_expr = pd.DataFrame({
    'GENE': df1['GENE'],
    't_statistic': t_statistics,
    'p_value': p_values
})

# Calculate the absolute value of the t-statistic for each gene
diff_expr['abs_t_statistic'] = np.abs(diff_expr['t_statistic'])

# Sort the genes by the absolute value of the t-statistic (in descending order)
diff_expr = diff_expr.sort_values(by='abs_t_statistic', ascending=False)

# Get the top 40 most differentially expressed genes
top_diff_expr = diff_expr.head(40)

top_diff_expr

# Write the differential expression results to .csv files
diff_expr.to_csv(f"{Clust_A}v{Clust_B}_diff_expr.csv", index=False)
top_diff_expr.to_csv(f"{Clust_A}v{Clust_B}_top_diff_expr.csv", index=False)


# Sort output by t-statistic and save list to .txt file to upload to single cell portal for viewing 
output_file_path = f"{Clust_A}v{Clust_B}_top40.txt"

sorted_top_diff_expr = top_diff_expr.sort_values(by=top_diff_expr.columns[1])
# Write the first column aka gene names to the text file
sorted_top_diff_expr.iloc[:, 0].to_csv(output_file_path, index=False, header=False)







