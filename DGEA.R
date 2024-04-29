# Load necessary libraries
library(scater)
library(scran)
library(igraph)
library(SingleCellExperiment)
library(fgsea)
library(msigdb)
library(ExperimentHub)
library(ggplot2)
library(reshape2)


# Set the working directory to the location of your data files.
setwd("/home/ajo26/SRP") 

# Read in your data
cc_logcounts <- read.csv("normalized_log_counts.csv", row.names = 1)

# Convert the gene expression data to a matrix
cc_logcounts_matrix <- as.matrix(cc_logcounts)

# Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(logcounts = cc_logcounts_matrix))

# Perform dimensionality reduction
sce <- runPCA(sce)

# Run t-SNE on the PCA-reduced data.
set.seed(1000) 
sce <- runTSNE(sce, dimred = "PCA", perplexity=50, n_dimred=10)

# Construct a shared nearest neighbor graph
snn_graph_tsne <- buildSNNGraph(sce, use.dimred = "TSNE")

# Perform clustering using the Louvain algorithm
clusters <- igraph::cluster_louvain(snn_graph_tsne)
colData(sce)$clusters <- factor(clusters$membership)

# Make the clusters variable  a factor 
clusters <- factor(colData(sce)$clusters)

#Visualize clustering of cell corrected data
plotReducedDim(sce, "TSNE", colour_by = "clusters")

# Identify marker genes using the findMarkers function
marker_genes <- findMarkers(sce, groups=clusters, test.type="t", 
                            pval.type="some", min.prop=0.2, 
                            log.p=TRUE, full.stats=FALSE)

# Create an ExperimentHub instance
eh <- ExperimentHub()

# Retrieve the MSigDB v7.2 human gene sets with gene symbols
msigdb.hs <- eh[["EH5422"]]

hallmarks <- subsetCollection(msigdb.hs, 'h')
c2 <- subsetCollection(msigdb.hs, 'c2')
c5 <- subsetCollection(msigdb.hs, 'c5')

# Combine the gene sets into one list
gene_sets <- c(hallmarks, c2, c5)

# Extract the summary.logFC for each cluster's markers
logFC_list <- lapply(marker_genes, function(df) {
  # Convert DFrame to a regular data.frame
  df <- as.data.frame(df)
  
  # Extract the summary.logFC vector
  fc_vector <- df$summary.logFC
  
  # Use gene names from rownames as names for the logFC vector
  names(fc_vector) <- rownames(df)
  
  return(fc_vector)
})

# Combine the gene sets into one list of character vectors
gene_sets <- c(geneIds(hallmarks), geneIds(c2), geneIds(c5))

# Prepare a list to store GSEA results
fgsea_results <- list()

# Perform GSEA on each cluster's logFC vector using fgseaMultilevel
for (i in seq_along(logFC_list)) {
  # Run fgseaMultilevel without the nperm argument
  results <- fgseaMultilevel(pathways = gene_sets, stats = logFC_list[[i]], minSize = 15, maxSize = 500)
  # Store results in the list
  fgsea_results[[i]] <- results
}

# Names the results list
names(fgsea_results) <- paste("Cluster", seq_along(logFC_list), sep = "_")

# Convert GSEA results to a data frame for plotting
gsea_df <- do.call(rbind, lapply(seq_along(fgsea_results), function(i) {
  res <- fgsea_results[[i]]
  if (is.null(res)) return(NULL)  # Skip if results are NULL
  res$Cluster <- paste("Cluster", i, sep = "_")
  return(res)
}))

# Filter for significant pathways if desired
gsea_df <- gsea_df[gsea_df$pval < 0.05 & !is.na(gsea_df$NES), ]

# Pivot the data frame into a wide format where rows are pathways and columns are clusters
gsea_matrix <- reshape2::dcast(gsea_df, pathway ~ Cluster, value.var = "NES")

# Set row names to pathways
rownames(gsea_matrix) <- gsea_matrix$pathway
gsea_matrix$pathway <- NULL  

# Remove pathways that have NA in any cluster from gsea_matrix before clustering
gsea_matrix <- na.omit(gsea_matrix)

# Ensure that gsea_matrix is a numeric matrix 
gsea_matrix_numeric <- as.matrix(gsea_matrix)
rownames(gsea_matrix_numeric) <- rownames(gsea_matrix)
gsea_matrix <- gsea_matrix_numeric

# Perform hierarchical clustering
row_dendrogram <- hclust(dist(gsea_matrix))
col_dendrogram <- hclust(dist(t(gsea_matrix)))

# Plot the heatmap
pheatmap::pheatmap(gsea_matrix, 
                   cluster_rows = row_dendrogram, 
                   cluster_cols = col_dendrogram, 
                   scale = "row",
                   color = colorRampPalette(c("blue", "white", "red"))(255), 
                   border_color = NA)  