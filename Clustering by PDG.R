if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("scater", "scran", "igraph", "Rtsne", "SingleCellExperiment"))

# Load necessary libraries 
library(SingleCellExperiment)
library(scran)
library(scater)
library(igraph)
library(ggplot2)

# Set the working directory to the location of your data files.
setwd("/home/ajo26/SRP") 

# Load the gene expression counts from a CSV file(exchange with cell cycle corrected counts after cell cycle correction)
counts <- read.csv("countsafterQC.csv", row.names = 1)

# Convert counts to matrix 
counts <- as.matrix(counts)

# Normalize each library to 100k reads. 
libSizes <- colSums(counts)
sizeFactors <- libSizes / 100000

# Apply the size factors to scale the counts.
counts_scaled <- t(t(counts) / sizeFactors)

# Add a pseudocount of one and log-transform.
logcounts <- log1p(counts_scaled)

# Create a SingleCellExperiment object.
sce <- SingleCellExperiment(assays = list(logcounts = logcounts))

# Calculate PDG for each cell (percentage of genes with detected expression).
pdg <- colSums(assay(sce, "logcounts") > 0) / nrow(sce) * 100
colData(sce)$PDG <- pdg

# Bin PDG into quartiles for visualization.
pdg_quartiles <- quantile(pdg, probs = seq(0, 1, 0.25))
pdg_bins <- cut(pdg, breaks = pdg_quartiles, labels = paste0("L", 1:4), include.lowest = TRUE)
colData(sce)$PDG_bins <- pdg_bins

# Perform PCA on the log-transformed normalized data.
sce <- runPCA(sce, exprs_values = "logcounts")

# Run t-SNE on the PCA-reduced data.
set.seed(1000) 
sce <- runTSNE(sce, dimred = "PCA", perplexity=50, n_dimred=10)

# Compute nearest-neighbour graphs using t-SNE and UMAP embeddings with scran.
snn_graph_tsne <- buildSNNGraph(sce, use.dimred = "TSNE")

# Perform graph-based clustering using the Louvain algorithm on the t-SNE graph.
cluster_tsne <- cluster_louvain(snn_graph_tsne)

# Add clustering results to SCE.
colData(sce)$Cluster_tsne <- factor(cluster_tsne$membership)


# Extract the t-SNE coordinates and cluster assignments for plotting
tsne_data <- data.frame(
  TSNE1 = reducedDim(sce, "TSNE")[, 1],
  TSNE2 = reducedDim(sce, "TSNE")[, 2],
  Cluster = colData(sce)$Cluster_tsne
)

# Using ggplot2 to plot the t-SNE colored by the clusters
p_clusters <- ggplot(tsne_data, aes(x = TSNE1, y = TSNE2, color = Cluster)) +
  geom_point(alpha = 0.5) +
  scale_color_discrete(name = "Clusters") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()  
  )

# Print the plot
print(p_clusters)

# Visualize t-SNE with cells colored by PDG.
p_pdg <- plotTSNE(sce, colour_by = "PDG_bins") + 
  labs(title = "Percentage of Detected Genes", colour = "PDG") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Print the plot
print(p_pdg)

# Prepare data for the bar plot
# Creating data frame that aggregates the count of each PDG bin by cluster
pdg_data <- data.frame(
  Cluster = colData(sce)$Cluster_tsne,
  PDG_Bins = colData(sce)$PDG_bins
)

# Calculate the counts of each PDG bin within each cluster
pdg_counts <- table(pdg_data)

# Turn into a proportion
pdg_props <- prop.table(pdg_counts, margin = 1)

# Melt the table into a long format for ggplot2
pdg_long <- as.data.frame(as.table(pdg_props))
colnames(pdg_long) <- c("Cluster", "PDG_Bins", "Proportion")

# Create the bar plot
pdg_bar_plot <- ggplot(pdg_long, aes(x = Cluster, y = Proportion, fill = PDG_Bins)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cluster", y = "Cluster Composition (%)", fill = "PDG") +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Print the bar plot
print(pdg_bar_plot)
