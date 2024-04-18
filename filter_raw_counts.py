#!/usr/bin/env python
# coding: utf-8

import pandas as pd


# Load the data
count_matrix = pd.read_csv('GSE202695_counts_raw.csv', index_col=0)
fastqscreen_counts = pd.read_csv('fastqscreen_counts.csv')
cell_metadata = pd.read_csv('cell_metadata.csv')


# Merge fastqscreen counts with cell metadata
metadata_with_counts = pd.merge(cell_metadata, fastqscreen_counts, on='Run_ID')

# Filter based on human-to-mouse ratio > 5 and human reads > 100,000
valid_metadata = metadata_with_counts[
    (metadata_with_counts['Human'] / metadata_with_counts['Mouse'] >= 5) &
    (metadata_with_counts['Human'] >= 100000)
]

# Get the list of valid Cell_IDs from the filtered metadata
valid_cell_ids = valid_metadata['Cell_ID'].tolist()

# Remove doublet and debris cells from valid Cell_IDs list
valid_cell_ids = [cell_id for cell_id in valid_cell_ids if '2cells' not in cell_id and 'debris' not in cell_id]


# Filter the count matrix to include only valid cells
filtered_count_matrix = count_matrix[valid_cell_ids]


# Compute the number of genes expressed per cell (non-zero counts)
genes_expressed_per_cell = (filtered_count_matrix > 0).sum(axis=0)

# Filter cells by gene count thresholds
filtered_cells_by_gene_count = genes_expressed_per_cell[
    (genes_expressed_per_cell >= 2346) & (genes_expressed_per_cell <= 9884)
].index

# Retain only columns (cells) that meet the gene expression criteria
final_filtered_count_matrix = filtered_count_matrix[filtered_cells_by_gene_count]

# Save the final filtered count matrix to a new CSV file
final_filtered_count_matrix.to_csv('countsafterQC.csv')






