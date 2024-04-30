library(Seurat)
library(SeuratObject)
library(fgsea)
library(msigdb)
library(ExperimentHub)
library(dplyr)
library(reshape2)
library(pheatmap)
library(GSEABase)

# Set the working directory to the location of your data files.
setwd("/home/ajo26/SRP") 

# Read and prepare data
raw_data <- read.csv("GSE202695_counts_afterQC.csv", row.names = 1, check.names = FALSE)
data_matrix <- as.matrix(raw_data)

# Create and preprocess Seurat object
seurat_object <- CreateSeuratObject(counts = data_matrix)
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object, features = VariableFeatures(object = seurat_object))

# Load the Seurat provided cell cycle gene sets
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Score cell cycle phases
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes)

# Run PCA as prerequisite for t-SNE
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Run t-SNE
seurat_object <- RunTSNE(seurat_object, dims = 1:10)

# Find neighbors and perform clustering
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Run t-SNE and plot with cell cycle phases before regressing out
seurat_object <- RunTSNE(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "tsne", group.by = "Phase", label = TRUE)

# Score cell cycle phases
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes)

# Calculate the difference between G2M and S phase scores
seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score

# Regress out the difference in cell cycle scores
seurat_object <- ScaleData(seurat_object, vars.to.regress = "CC.Difference", features = rownames(seurat_object))

# Compute PCA and t-SNE with the cell cycle differences regressed out
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- RunTSNE(seurat_object, dims = 1:10)

# Find neighbors and perform clustering with the updated data
seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# plotting after regressing
cell_cycle_genes <- c(s.genes, g2m.genes)
seurat_object <- RunPCA(seurat_object, features = cell_cycle_genes)
DimPlot(seurat_object, reduction = "pca", group.by = "Phase")

# After clustering
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# Access the cluster IDs
cluster_ids <- Idents(seurat_object)

# Identify differentially expressed genes for each cluster
all_markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)

# Inspect the column names of the all_markers data frame
colnames(all_markers)

# Filter for genes that are differentially expressed in at least 20% of comparisons
# Calculate the minimum p-value for each gene across comparisons
# Ensure that avg_log2FC is also included in the summarization
top_markers <- all_markers %>%
  group_by(gene) %>%
  summarize(
    min_p_val_adj = min(p_val_adj),
    avg_log2FC = mean(avg_log2FC) # Correct column name used here
  ) %>%
  # Select the top genes based on the min_p_val_adj
  top_n(n = n()/5, wt = min_p_val_adj)

# Arrange the top markers by average log fold change in descending order
ranked_genes <- top_markers %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

# Initialize ExperimentHub to access MSigDB datasets
eh <- ExperimentHub()

# Retrieve the MSigDB v7.2 human gene sets with gene symbols
msigdb.hs <- eh[["EH5421"]]

# Convert GeneSetCollection to a list of gene sets
# This extracts each GeneSet, converts it to a character vector of genes
gene_sets_list <- lapply(as.list(msigdb.hs), function(gs) {
  as.character(geneIds(gs))
})

# Load collections
hallmarks <- subsetCollection(msigdb.hs, 'h')
c2 <- subsetCollection(msigdb.hs, 'c2')
c5 <- subsetCollection(msigdb.hs, 'c5')

# Combine the gene sets into one list correctly
hallmarks_list <- lapply(as.list(hallmarks), function(gs) as.character(geneIds(gs)))
c2_list <- lapply(as.list(c2), function(gs) as.character(geneIds(gs)))
c5_list <- lapply(as.list(c5), function(gs) as.character(geneIds(gs)))

# Combine all gene sets into one list
all_gene_sets <- c(hallmarks_list, c2_list, c5_list)

# Arrange the top markers by average log fold change in descending order
ranked_genes <- top_markers %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)

# Assign the 'pathway' column to each gene set
all_gene_sets <- lapply(all_gene_sets, function(set) {
  if (length(set) > 0) {
    data.frame(pathway = rep(names(set), each = length(set)), genes = unlist(set))
  } else {
    data.frame(pathway = character(), genes = character())
  }
})

# Create the named stats vector with gene names as names and avg_log2FC as the values
stats <- setNames(as.numeric(top_markers$avg_log2FC), top_markers$gene)


# Run fgsea
fgsea_results <- fgsea(pathways = all_gene_sets,
                       stats = stats,
                       minSize = 15,
                       maxSize = 500,
                       scoreType = "pos") 


