#### This R script is designed to replicate Figures 1D and 1E from the cited article, focusing on the heterogeneity 
# in metastatic breast cancer through single-cell RNA sequencing analysis. The analysis workflow includes:
#   - Performing PCA to reduce the dimensionality of the data, capturing the main variations across cells.
#   - Applying t-SNE to visualize these variations in a two-dimensional plot, aiding in the identification of cellular clusters.
#   - Constructing an SNN graph from PCA results to define neighbourhoods based on shared nearest neighbours.
#   - Using the Louvain method to cluster cells within the SNN graph, aiming to detect clusters that correspond to different cell types or conditions.
#   - Visualizing the identified clusters and comparing them with known cell types or lines using t-SNE plots.

#### Purpose:
# The script demonstrates an attempt to independently replicate the clustering results presented in the original study by 
# applying the same data preprocessing and analysis steps. This ensures the reproducibility of the results and validates the 
# methodology described by the original authors.



# Load necessary libraries for analysis and visualization.
library(SingleCellExperiment)
library(scran)
library(scater)
library(igraph)
library(ggplot2)

# Define the working directory containing the data files.
setwd("/alice-home/1/e/es452/Research_Project/AleroData/")

# Import the count matrix and cell metadata, setting the first column as row names.
data <- read.csv("countsafterQC.csv", row.names = 1)
metadata <- read.csv("filtered_cell_metadata.csv", row.names = 1)

# Construct a SingleCellExperiment object with the imported count matrix.
sce <- SingleCellExperiment(assays = list(counts = as.matrix(data)))

# Normalize the count data using log-transformation.
sce <- logNormCounts(sce)

# Carry out Principal Component Analysis (PCA) to reduce data dimensionality.
sce <- runPCA(sce, ncomponents = 50)

# Perform t-SNE to further reduce dimensions and allow for visualization.
sce <- runTSNE(sce, dimred = "PCA")

# Construct a Shared Nearest Neighbor (SNN) graph from the PCA results.
snn_graph <- buildSNNGraph(sce, k = 5, use.dimred = 'PCA')

# Use the Louvain algorithm to detect communities/clusters in the SNN graph.
clusters <- igraph::cluster_louvain(snn_graph)

# Convert the 'cell_line' column to a factor for categorical plotting.
metadata$cell_line <- factor(metadata$cell_line)

# Append the clustering results to the metadata for later use in plotting.
metadata$Cluster <- factor(clusters$membership)

# Extract t-SNE coordinates from the SingleCellExperiment object.
tsne_coords <- reducedDims(sce)$TSNE
tsne_data <- as.data.frame(tsne_coords)

# Align 'cell_line' data from metadata with the t-SNE data.
tsne_data$cell_line <- metadata$cell_line[match(rownames(tsne_data), rownames(metadata))]

# Remove rows with missing 'cell_line' information to ensure clean plotting.
tsne_data <- tsne_data[!is.na(tsne_data$cell_line), ]

### Creating Figure 1d: t-SNE Plot of Sequenced Cells by Cell Line.

# Confirm that 'cell_line' is treated as a categorical variable in plots.
tsne_data$cell_line <- factor(tsne_data$cell_line)

# Define a mapping from cell line codes to more descriptive labels.
name_mapping <- c("HBRX1921" = "PDX1",
                  "HBRX2344" = "PDX4",
                  "HBRX2353" = "PDX2",
                  "HBRX3078" = "PDX3",
                  "MDAMB231" = "MDAMB231")

colors <- c("PDX1" = "green", "PDX2" = "purple", "PDX3" = "red", 
            "PDX4" = "blue", "MDAMB231" = "orange")

# Apply the mapping to the 'cell_line' column for clearer visualization.
tsne_data$cell_line <- factor(name_mapping[tsne_data$cell_line])

# Generate the t-SNE plot, colored by 'cell_line', and display it.
tsne_plot_by_cell_line <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = cell_line)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = colors) +
  labs(title = "t-SNE Plot of Sequenced Cells by Cell Line",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic() +
  theme(legend.position = "right")
print(tsne_plot_by_cell_line)


### Creating Figure 1e: t-SNE Plot of Sequenced Cells by Type (Tissue)

# Prepare 'tissue' column for plotting, similar to 'cell_line'.
metadata$tissue <- factor(metadata$tissue)
tsne_data$tissue <- metadata$tissue[match(rownames(tsne_data), rownames(metadata))]
tsne_data <- tsne_data[!is.na(tsne_data$tissue), ]
tsne_data$tissue <- factor(tsne_data$tissue)

# Generate the t-SNE plot, colored by 'tissue', and display it.
tsne_plot_by_type <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = tissue)) +
  geom_point(size = 1, alpha = 0.8) +  # Adjust size and transparency as needed
  scale_color_manual(values = c("Metastases" = "blue", "Tumor" = "orange")) +  # Use appropriate colors
  labs(title = "t-SNE Plot of Sequenced Cells by Origin",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic() 
print(tsne_plot_by_type)
