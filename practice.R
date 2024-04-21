library(Seurat)
library(ggplot2)
library(dplyr)
# Load the data
data <- read.csv("~/practice_sre/GSE202695_counts_afterQC.csv", row.names = 1)
metadata <- read.csv("~/practice_sre/GSE202695_metadata.csv", row.names = 1)
# Create Seurat object with data and metadata
seurat_obj <- CreateSeuratObject(counts = data)
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
# Normalize data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
# Perform PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
#Run t-SNE using PCA results
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)
# Build SNN graph based on the PCA results
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
# Perform Louvain clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
# Extract t-SNE coordinates
# Access t-SNE coordinates stored in the 'reductions' slot
tsne_coords <- Embeddings(seurat_obj, reduction = "tsne")
# Convert the matrix to a dataframe
tsne_data <- as.data.frame(tsne_coords)
# Rename the columns appropriately
colnames(tsne_data) <- c("TSNE1", "TSNE2")
# Add metadata to the tsne_data dataframe
tsne_data$cluster <- Idents(seurat_obj)  # This adds the cluster information
# Add cluster and model information to the t-SNE data frame
tsne_data$Cluster <- Idents(seurat_obj)
tsne_data$model <- seurat_obj$model
# Map the 'model' column to new names using the mapping vector
name_mapping <- c("HBRX1921" = "PDX1", "HBRX2344" = "PDX4", "HBRX2353" = "PDX2", "HBRX3078" = "PDX3", "MDAMB231" = "MDAMB231")
tsne_data$model <- factor(name_mapping[tsne_data$model])
# Assuming you have these models, specify colors for each.
colors <- c("PDX1" = "green", "PDX2" = "purple", "PDX3" = "red", "PDX4" = "blue", "MDAMB231" = "orange")
# Plot the t-SNE
tsne_plot <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = tsne_data$model)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = colors) +
  labs(title = "t-SNE Plot of Sequenced Cells by Model", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_classic()
print(tsne_plot)
# Assuming 'type' is in your metadata
tsne_data$type <- factor(seurat_obj$type)
# Check the column names of tsne_data
colnames(tsne_data)
tsne_plot_by_type <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = type)) +
  geom_point(size = 1.5, shape = 16) +  # Customize point size and shape
  scale_color_brewer(palette = "Set1") +  # Use a color palette that improves visibility
  labs(title = "t-SNE Plot of Sequenced Cells by Type",
       x = "t-SNE Dimension 1", y = "t-SNE Dimension 2") +
  theme_minimal()  # A minimal theme for cleaner look


