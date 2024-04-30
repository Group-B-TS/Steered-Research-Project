# Required libraries
library(devtools)
library(griph)         # For cell cycle prediction
library(org.Hs.eg.db)  # Gene annotations
library(AnnotationDbi) # Annotation conversion
library(dplyr)         # Data manipulation

# Install and load these if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")

# Load data
raw_data <- read.csv("/home/im205/scRNAseq_project/countsafterQC.csv", row.names = 1, check.names = FALSE)
data_matrix <- as.matrix(raw_data)

# Convert gene symbols to Entrez IDs
entrez_ids <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = rownames(data_matrix),
                                    column = "ENTREZID",
                                    keytype = "SYMBOL",
                                    multiVals = "first")

# Filter valid genes
valid_ids <- !is.na(entrez_ids)
data_matrix <- data_matrix[valid_ids, ]
rownames(data_matrix) <- entrez_ids[valid_ids]

cell_cycle_scores <- griph::predictCellCycle(cnt = data_matrix, org = "human.Whitfield", cor_thr = 0.2, granularity = "high", refine_iter = 200)

# Assuming cell_cycle_scores are correctly computed and are aligned with the columns of data_matrix
if (length(cell_cycle_scores) != ncol(data_matrix)) {
    stop("Mismatch in the number of cells and cell cycle scores")
}

# Create a data frame with necessary variables for GLM
df_for_glm <- data.frame(
    cell_cycle = factor(cell_cycle_scores),  # Convert cell cycle scores to a factor if they are categorical
    total_counts = log(colSums(data_matrix)) # Library size normalization factor
)

# Generate the model matrix for GLM, excluding the intercept if needed
model_matrix <- model.matrix(~ cell_cycle + total_counts - 1, data = df_for_glm)


# Check the model matrix
head(model_matrix)

# Assuming data_matrix is normalized
corrected_expression <- matrix(nrow = nrow(data_matrix), ncol = ncol(data_matrix), dimnames = dimnames(data_matrix))

# Apply GLM to each gene's expression
for (i in 1:nrow(data_matrix)) {
    fit <- glm(data_matrix[i, ] ~ cell_cycle + total_counts - 1, data = df_for_glm, family = gaussian())
    corrected_expression[i, ] <- residuals(fit)
}
corrected_expression[corrected_expression <= 0] <- 1e-4  # Replace non-positive values to a small positive number before log transformation

# Check results
dim(corrected_expression)  # Should match dimensions of the original data_matrix

# Add an offset to ensure there are no negative values, if not already done
min_value <- min(corrected_expression)
if (min_value <= 0) {
    corrected_expression <- corrected_expression + abs(min_value) + 1
}

# Normalize expression data
normalized_expression <- log1p(corrected_expression)  # log1p is safe for non-negative values

write.csv(normalized_expression, "/home/im205/scRNAseq_project/normalized_log_counts.csv", row.names = TRUE)

# Assuming 'normalized_expression' is your data matrix
# Convert to a matrix if not already (ensure it's in the correct format)
normalized_matrix <- as.matrix(normalized_expression)

# Remove duplicate rows
unique_normalized_matrix <- normalized_expression[!duplicated(normalized_expression), ]
# Check how many duplicates were removed
cat("Removed", nrow(normalized_matrix) - nrow(unique_normalized_matrix), "duplicate rows\n")

# Now, perform t-SNE on the unique data
tsne_result <- Rtsne(unique_normalized_matrix, dims = 2, perplexity = 30, verbose = TRUE)
tsne_data <- data.frame(tsne_result$Y)

# Calculate variance for each gene
gene_variances <- apply(normalized_expression, 1, var)

# Filter out genes with zero variance
filtered_expression <- normalized_expression[gene_variances > 0, ]

# Confirm that filtered_expression is correctly calculated
if (ncol(filtered_expression) == 0) {
    stop("No variable genes found after filtering.")
}
# Check if there are genes left after filtering
if (nrow(filtered_expression) > 0) {
    pca_result <- prcomp(t(filtered_expression), scale. = TRUE)
    
    # Plot PCA results
    pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2])
    library(ggplot2)
    ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point() +
        theme_minimal() +
        ggtitle("PCA of Filtered Expression Data")
} else {
    stop("No variable genes available for PCA.")
}

# Hierarchical clustering
dist_matrix <- dist(t(unique_normalized_matrix))
hc <- hclust(dist_matrix)  # Hierarchical clustering

# Cut the dendrogram to create clusters
clusters <- cutree(hc, k = 5)  # Adjust k based on your data

tsne_data$cluster <- as.factor(clusters[!duplicated(normalized_expression)])

library_complexity <- log(colSums(exp(normalized_expression) - 1))

# Plotting
tsne_result <- Rtsne(unique_normalized_matrix, dims = 2, perplexity = 30, verbose = TRUE)
tsne_data <- data.frame(tsne_result$Y)
library_complexity <- log(colSums(data_matrix))  # Assure this matches the filtered/corrected data matrix
data <- data.frame(
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2],
    cluster = tsne_data$cluster,
    library_complexity = library_complexity,  # Match length with tsne_result$Y
    cell_cycle_stage = factor(cell_cycle_scores[!duplicated(tsne_result$Y)])  # Ensure correct indexing
)
data <- data[!is.na(data$cell_cycle_stage), ]
data$library_complexity_log <- log1p(data$library_complexity)

#2a Library Complexity
ggplot(data, aes(x = tSNE1, y = tSNE2, color = library_complexity_log)) +
    geom_point(alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "orange") +  # Adjust colors as needed
    labs(title = "Library Complexity (log-transformed)", x = "t-SNE 1", y = "t-SNE 2") +
    theme_minimal() +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))

# 2b Cell Cycle Similarity
ggplot(data, aes(x = tSNE1, y = tSNE2, color = cell_cycle_stage)) +
    geom_point(alpha = 0.6) +
    labs(title = "Cell Cycle Similarity") +
    theme_minimal()

#2c t-SNE Clustering of Cells
ggplot(data, aes(x = tSNE1, y = tSNE2, color = cluster)) +
    geom_point(alpha = 0.6) +
    labs(title = "t-SNE Clustering of Cells") +
    theme_minimal() +
    scale_color_viridis_d()

