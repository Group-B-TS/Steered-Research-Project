library(Seurat)
# Load data
counts <- read.csv(file = "~/practice_sre/GSE202695_counts_afterQC.csv")
metadata <- read.csv("~/practice_sre/transposed_cell_data.csv", header = TRUE, row.names = 1)
# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj
Assays(seurat_obj)
seurat_obj
seurat_obj@assays$RNA$counts[is.na(seurat_obj@assays$RNA$counts)] <- 0
# Normalization and scaling
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj
# Further filtering based on gene expression
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 2346 & nFeature_RNA < 9884)
# Calculate the total number of unique genes
total_genes <- length(row.names(seurat_obj))
# Calculate thresholds for 10% and 40%
lower_threshold <- total_genes * 0.10
upper_threshold <- total_genes * 0.40
# Filter cells based on these thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= lower_threshold & nFeature_RNA <= upper_threshold)
# Check the results
print(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
# Run t-SNE using the results of PCA
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)
# Add metadata to the Seurat object
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
# Define color mapping for plotting
name_mapping <- c("HBRX1921" = "PDX1", "HBRX2344" = "PDX4", "HBRX2353" = "PDX2", 
                  "HBRX3078" = "PDX3", "MDAMB231" = "MDAMB231")
colors <- c("PDX1" = "green", "PDX2" = "purple", "PDX3" = "red", "PDX4" = "blue", "MDAMB231" = "orange")

# Map metadata 'model' to new names
seurat_obj$cell_line <- name_mapping[seurat_obj$cell_line]

# Plot t-SNE
DimPlot(seurat_obj, reduction = "tsne", group.by = "cell_line", label = TRUE, label.size = 3, pt.size = 1) +
  scale_color_manual(values = colors) +
  labs(title = "t-SNE Plot of Sequenced Cells by Model of Origin") +
  theme_classic()
# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
# Plot UMAP
DimPlot(seurat_obj, reduction = "umap")

