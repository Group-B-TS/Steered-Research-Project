###                   Overview:             

### This R script is designed to replicate Figures 1D and 1E from the cited article, focusing on the heterogeneity 
# in metastatic breast cancer through single-cell RNA sequencing analysis.
# The script demonstrates the robustness of the Seurat package in performing comprehensive single-cell analysis, with a focus on:
#   - Using the exact filtered single-cell RNA-seq data as provided by the authors of the study.
#   - Normalizing the dataset to facilitate accurate comparisons across single cells.
#   - Executing PCA to distill the dataset into its principal components, thus simplifying the complexity inherent to single-cell data.
#   - Employing t-SNE as a method for dimensionality reduction, aiming to visualize cellular clusters in a two-dimensional space.
#   - Visualizing the cell clusters in t-SNE plots, enriched with annotations to reflect known cell types or lines.

#### Data Source: 
# The script processes the exact same filtered dataset made publicly available by the original authors, which can be retrieved from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202695

### Purpose:
# The primary goal of this script is to challenge and verify the reproducibility of published findings by employing novel bioinformatics pipelines.
# It is crafted to assess whether different computational strategies, such as those embedded in the Seurat package, lead to similar interpretations 
#of cellular heterogeneity as those reported in the original study.

### Script Origin:
# This code was structured upon and inspired by the "Zainab Practice.R" and "Eleftherios t-SNE Clusterig.R" script.

# Load necessary libraries for Seurat analysis and data visualization
library(Seurat)
library(ggplot2)
library(dplyr)

# Set the working directory to the location of the data files
setwd("/alice-home/1/e/es452/Research_Project/Zainab/")

# Load the gene expression counts and metadata into R
data <- read.csv("GSE202695_counts_afterQC.csv", row.names = 1)
metadata <- read.csv("GSE202695_metadata.csv", row.names = 1)

# Create a Seurat object with the raw count data
seurat_obj <- CreateSeuratObject(counts = data)

# Incorporate metadata into the Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

# Normalize the data using a log-normalization method
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify variable features likely to be informative
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data to have zero mean and unit variance, preparing for PCA
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# Run PCA using the variable features to reduce dimensionality
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Execute t-SNE using the results of PCA (first 20 dimensions)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)

# Find neighbors to construct a graph for clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Cluster cells using the Louvain algorithm based on the graph constructed above
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Extract t-SNE coordinates from the Seurat object
tsne_coords <- Embeddings(seurat_obj, reduction = "tsne")

# Convert the matrix to a data frame for use with ggplot2
tsne_data <- as.data.frame(tsne_coords)
colnames(tsne_data) <- c("TSNE1", "TSNE2")

### Creating Figure 1d: t-SNE Plot of Sequenced Cells by Cell Line ###

# Annotate t-SNE data with cluster assignments and model information from metadata
# Ensuring that cluster and model information are aligned with the t-SNE data
tsne_data$Cluster <- Idents(seurat_obj)
tsne_data$model <- metadata$model[match(rownames(tsne_data), rownames(metadata))]

# Filter out NAs in the model column before plotting
tsne_data <- tsne_data[!is.na(tsne_data$model), ]

# Create a name mapping from raw to more interpretable model names
name_mapping <- c("HBRX1921" = "PDX1", "HBRX2344" = "PDX4", "HBRX2353" = "PDX2", 
                  "HBRX3078" = "PDX3", "MDAMB231" = "MDAMB231")

# Apply the name mapping to the 'model' column
tsne_data$model <- factor(name_mapping[tsne_data$model])

# Define a color palette for the different models and plot the t-SNE with colors mapped to models
colors <- c("PDX1" = "green", "PDX2" = "purple", "PDX3" = "red", 
            "PDX4" = "blue", "MDAMB231" = "orange")

# Plot the t-SNE, now with no NAs and with a color legend for each model
tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = model)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = colors) +
  labs(title = "t-SNE Plot of Sequenced Cells by Model", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic() +
print(tsne_plot)

### Creating Figure 1e: t-SNE Plot of Sequenced Cells by Type ###

# Assuming 'type' is present in your metadata and you want to plot by type as well
tsne_data$type <- metadata$type[match(rownames(tsne_data), rownames(metadata))]

# Filter out NAs in the type column to prevent them from appearing in the plot
tsne_data <- tsne_data[!is.na(tsne_data$type), ]
tsne_data$type <- factor(tsne_data$type)

# Plot the t-SNE color-coded by cell type, ensuring that type information is correctly mapped
tsne_plot_by_type <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = type)) +
  geom_point(size = 1.5, shape = 16) +
  scale_color_manual(values = c("Lung" = "blue", "Tumor" = "orange")) +
  labs(title = "t-SNE Plot of Sequenced Cells by Type", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic() +
print(tsne_plot_by_type)
