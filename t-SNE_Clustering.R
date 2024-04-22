### This R script is designed to replicate Figures 1D and 1E from the cited article, focusing on the heterogeneity 
# in metastatic breast cancer through single-cell RNA sequencing analysis. The analysis workflow includes:
#   - Utilizing the exact filtered single-cell RNA-seq data provided by the authors.
#   - Normalizing the data and performing PCA to reduce the dimensionality of the data.
#   - Applying t-SNE to project these variations into a two-dimensional plot for visualization and to aid in identifying cellular clusters.
#   - Constructing an SNN graph from PCA results to define neighborhoods based on shared nearest neighbors.
#   - Using the Louvain method to detect communities within the SNN graph, which correspond to different cell types or conditions.
#   - Visualizing the identified clusters and comparing them with known cell types or lines using t-SNE plots to observe patterns of heterogeneity.

#### Data source: Data used in this script is the exact same filtered data as provided by the authors, available at
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202695

### Purpose:
# The purpose of this script is to verify the reproducibility of published results by using the authors' processed data and 
# following their analytical methodology as described. This provides a validation of both the data and the techniques used in the
# original study, ensuring that the conclusions drawn are robust and reliable.



# Load necessary libraries for handling single-cell data and visualization.
library(SingleCellExperiment)
library(scran)
library(scater)
library(igraph)
library(ggplot2)

# Specify the directory where data files are stored. This sets the path so all file operations are relative to this directory.
setwd("/alice-home/1/e/es452/Research_Project/")

# Load the gene expression counts and the corresponding metadata from CSV files.
data <- read.csv("GSE202695_counts_afterQC.csv", row.names = 1)
metadata <- read.csv("GSE202695_metadata.csv", row.names = 1)

# Create a SingleCellExperiment object to store and manage the gene expression data efficiently.
sce <- SingleCellExperiment(assays = list(counts = as.matrix(data)))

# Normalize the gene expression data to make it suitable for downstream analysis.
sce <- logNormCounts(sce)

# Perform Principal Component Analysis (PCA) to reduce the dimensionality of the data while preserving variance.
sce <- runPCA(sce, ncomponents = 50)

# Apply t-SNE (t-Distributed Stochastic Neighbor Embedding) to project the PCA-reduced data into a two-dimensional space for visualization.
sce <- runTSNE(sce, dimred = "PCA")

# Build a Shared Nearest Neighbor (SNN) graph from the PCA results. This graph helps to identify local structures within the data.
snn_graph <- buildSNNGraph(sce, k = 5, use.dimred = 'PCA')

# Cluster the cells using the Louvain algorithm, which detects communities in the SNN graph to identify densely connected groups of cells.
clusters <- igraph::cluster_louvain(snn_graph)

# Convert the 'model' column in the metadata to a factor for categorical analysis.
metadata$model <- factor(metadata$model)

# Append the clustering results as a new 'Cluster' column in the metadata.
metadata$Cluster <- factor(clusters$membership)

# Extract the two-dimensional t-SNE coordinates for plotting.
tsne_coords <- reducedDims(sce)$TSNE
tsne_data <- as.data.frame(tsne_coords)
tsne_data$Cluster <- metadata$Cluster[match(rownames(tsne_data), rownames(metadata))]


### Creating Figure 1d: t-SNE Plot of Sequenced Cells by Cell Line ###

# Match the 'model' information from metadata to the t-SNE data, ensuring each point in the plot corresponds to a specific model.
tsne_data$model <- metadata$model[match(rownames(tsne_data), rownames(metadata))]
tsne_data <- tsne_data[!is.na(tsne_data$model), ]

# Define a mapping from model identifiers to more descriptive names for better visualization clarity.
name_mapping <- c("HBRX1921" = "PDX1",
                  "HBRX2344" = "PDX4",
                  "HBRX2353" = "PDX2",
                  "HBRX3078" = "PDX3",
                  "MDAMB231" = "MDAMB231")

colors <- c("PDX1" = "green", "PDX2" = "purple", "PDX3" = "red", 
           "PDX4" = "blue", "MDAMB231" = "orange")

# Apply the descriptive name mapping to the 'model' column in the plot data.
tsne_data$model <- factor(name_mapping[tsne_data$model])

# Generate and display a t-SNE plot colored by 'model', illustrating the distribution of different cell models in the projected space.
tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = model)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = colors) +
  labs(title = "t-SNE Plot of Sequenced Cells by Model",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic()
print(tsne_plot)


### Creating Figure 1e: t-SNE Plot of Sequenced Cells by Type ###

# Align 'Type' information from metadata with the t-SNE plot data for visualization.
tsne_data$Type <- metadata$type[match(rownames(tsne_data), rownames(metadata))]

# Ensure 'Type' is treated as a categorical variable for appropriate coloring in the plot.
tsne_data$Type <- factor(tsne_data$Type)

# Generate and display a second t-SNE plot, now colored by 'Type', to show differences between cell types like Lung and Tumor.
tsne_plot_by_type <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = Type)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("Lung" = "blue", "Tumor" = "orange")) +
  labs(title = "t-SNE Plot of Sequenced Cells by Type",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic() 
print(tsne_plot_by_type)

  theme_classic() 
print(tsne_plot_by_type)
