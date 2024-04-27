library(Seurat)
library(SeuratObject)
library(org.Hs.eg.db) # Assuming this is installed via Bioconductor
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(griph)  # Ensure griph is installed and loaded

# Install Bioconductor if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

raw_data <- read.csv("/home/im205/scRNAseq_project/GSE202695_counts_afterQC.csv", row.names = 1, check.names = FALSE)
data_matrix <- as.matrix(raw_data)

# Extract gene symbols from the data matrix
symbols <- rownames(data_matrix)

# Convert gene symbols to Entrez IDs
entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = symbols,
                                    column = "ENTREZID",
                                    keytype = "SYMBOL",
                                    multiVals = "first")

# Remove genes without valid Entrez IDs
valid_ids <- !is.na(entrez_ids)
data_matrix <- data_matrix[valid_ids, ]
rownames(data_matrix) <- entrez_ids[valid_ids]

# Create Seurat object
seurat_object <- CreateSeuratObject(counts = data_matrix)

# Normalize the data
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object, features = VariableFeatures(object = seurat_object))

# Run PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Run UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# Predict cell cycle stages using griph
cell_cycle_scores <- predictCellCycle(cnt = data_matrix, org = "human.Whitfield", cor_thr = 0.2, granularity = "high", refine_iter = 200)

# Add cell cycle information to Seurat object metadata
seurat_object$cell_cycle <- as.factor(cell_cycle_scores)

# Regress out cell cycle effects from the data
seurat_object <- ScaleData(seurat_object, vars.to.regress = "cell_cycle")

# Visualize UMAP with cell cycle stages
DimPlot(seurat_object, reduction = "umap", group.by = "cell_cycle", label = TRUE)

# Extract cell identifiers with corresponding cell cycle stages from meta.data
if ("cell_cycle" %in% names(seurat_object@meta.data)) {
    cell_ids <- rownames(seurat_object@meta.data)  # These are cell identifiers
    cell_cycle_stages <- seurat_object@meta.data$cell_cycle  # Extract cell cycle stages
    
    # Filter valid entries (non-NA)
    valid_entries <- !is.na(cell_cycle_stages)
    if (sum(valid_entries) > 0) {  # Check if there are any valid entries
        cell_cycle_data <- data.frame(CellID = cell_ids[valid_entries], CellCycleStage = cell_cycle_stages[valid_entries], stringsAsFactors = FALSE)
        
        # Write to CSV
        write.csv(cell_cycle_data, file = "Cell_CycleStages.csv", row.names = FALSE)
        print("CSV file has been successfully written with data.")
    } else {
        print("No valid entries to write to CSV.")
    }
} else {
    print("No cell cycle stage data found in Seurat object metadata.")
}

