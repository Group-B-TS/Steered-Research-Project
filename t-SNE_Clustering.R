library(SingleCellExperiment)
library(scran)
library(scater)
library(igraph)
library(ggplot2)

# Set the working directory to where the files are located
setwd("/alice-home/1/e/es452/Research_Project/")

# Read in the data
data <- read.csv("GSE202695_counts_afterQC.csv", row.names = 1)
metadata <- read.csv("GSE202695_metadata.csv", row.names = 1)

# Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(data)))

# Normalize counts
sce <- logNormCounts(sce)

# Perform PCA
sce <- runPCA(sce, ncomponents = 50)

# Run t-SNE on the PCA results
sce <- runTSNE(sce, dimred = "PCA")

# Build SNN graph based on the PCA results
snn_graph <- buildSNNGraph(sce, k = 5, use.dimred = 'PCA')

# Perform Louvain clustering
clusters <- igraph::cluster_louvain(snn_graph)

metadata$model <- factor(metadata$model)

# Add clustering results to metadata
metadata$Cluster <- factor(clusters$membership)

# Extract t-SNE coordinates for plotting
tsne_coords <- reducedDims(sce)$TSNE
tsne_data <- as.data.frame(tsne_coords)
tsne_data$Cluster <- metadata$Cluster[match(rownames(tsne_data), rownames(metadata))]

# Match the 'model' information in metadata with the corresponding rows in 'tsne_data'
tsne_data$model <- metadata$model[match(rownames(tsne_data), rownames(metadata))]
tsne_data <- tsne_data[!is.na(tsne_data$model), ]

name_mapping <- c("HBRX1921" = "PDX1",
                  "HBRX2344" = "PDX4",
                  "HBRX2353" = "PDX2",
                  "HBRX3078" = "PDX3",
                  "MDAMB231" = "MDAMB231")

### Creating Figure 1d: t-SNE Plot of Sequenced Cells by Model
# Map the 'model' column in tsne_data to the new names using the mapping vector
tsne_data$model <- factor(name_mapping[tsne_data$model])

# Plot the t-SNE with the updated model names
tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = model)) +
  geom_point(alpha = 0.8) +  # Adjust alpha for point transparency
  labs(title = "t-SNE Plot of Sequenced Cells by Model",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic()
print(tsne_plot)

### Creating Figure 1e: t-SNE Plot of Sequenced Cells by Type

# Map the 'Type' information from 'metadata' to 'tsne_data' using the rownames to align them
tsne_data$Type <- metadata$type[match(rownames(tsne_data), rownames(metadata))]

# Make sure that 'Type' is a factor for proper coloring
tsne_data$Type <- factor(tsne_data$Type)

# Now plot using the 'Type' column to color the points
tsne_plot_by_type <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = Type)) +
  geom_point(size = 1, alpha = 0.8) +  # Adjust size and transparency as needed
  scale_color_manual(values = c("Lung" = "blue", "Tumor" = "orange")) +  # Use appropriate colors
  labs(title = "t-SNE Plot of Sequenced Cells by Type",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic() 
print(tsne_plot_by_type)

