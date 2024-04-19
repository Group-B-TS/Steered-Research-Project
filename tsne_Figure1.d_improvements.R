# Based on Eleftherios' first commit for the script entitled 'RFigure1.d_clustering.R'.
# In this script I added steps of normalization and some preprocessing as well as 
# experimented with perplexity values and learning rates compared to his first commit.
# This script is mostly a draft and for testing, many values and names need adjusting!
# This script was used to produce a final and refined version of Figure 1.d.

# Set working directory
setwd("/home/ib156/Desktop/")

# Loading necessary libraries
library(Rtsne)
library(ggplot2)
library(caret)

# Reading data
data <- read.csv("GSE202695_counts_afterQC.csv", header = TRUE)
meta <- read.csv("GSE202695_metadata.csv", header = TRUE)

# Setting gene names as row names and remove the first column from the data
rownames(data) <- data[, 1]
data <- data[,-1]

# Normalizing each column (sample) to sum to 100,000 reads and log-transforming with pseudocount of 1
normalized_data <- sweep(data, 2, colSums(data), FUN = "/") * 100000
log_transformed_data <- log2(normalized_data + 1)

# Transposing the log-transformed data so that genes are columns and samples are rows
data_t <- t(log_transformed_data)

# Removing columns with zero or near-zero variance
nzv <- nearZeroVar(data_t)
data_t <- data_t[, -nzv]

# Performing PCA for dimensionality reduction on the cleaned data
pca_results <- prcomp(data_t, scale. = TRUE)

# Selecting the number of principal components to keep (50 in this case)
pca_data <- as.data.frame(pca_results$x[,1:50])

# Experimenting with different perplexity values and learning rates
perplexities <- c(5, 30, 50)  # To adjust, this is just to test
learning_rates <- c(100, 200, 500)  # To adjust, this is just to test

# Function to run t-SNE
run_tsne <- function(pca_data, perplexity, learning_rate) {
  set.seed(40)  # Setting a seed for reproducibility but likely needs adjustment
  Rtsne(pca_data, dims = 2, perplexity = perplexity, 
        learning_rate = learning_rate, verbose = TRUE, max_iter = 1000)
}

# Running t-SNE with different parameters and storing results
tsne_results_list <- list()
for (p in perplexities) {
  for (lr in learning_rates) {
    tsne_key <- paste("perp", p, "lr", lr, sep = "_")
    tsne_results_list[[tsne_key]] <- run_tsne(pca_data, p, lr)
  }
}

# We should choose the t-SNE result that visually matches best... testing as below for now
chosen_tsne_key <- "perp_30_lr_200"
chosen_tsne <- tsne_results_list[[chosen_tsne_key]]

# Combining t-SNE output with the metadata
tsne_data <- as.data.frame(chosen_tsne$Y)
tsne_data$Model <- meta[, 'model']  # Name needs to be adjusted to correct column... 

# Plotting the t-SNE results
tsne_plot <- ggplot(tsne_data, aes(x = V1, y = V2, color = Model)) +
  geom_point() +
  labs(title = "t-SNE Plot of Sequenced Cells by Model of Origin",
       x = "t-SNE Dimension 1",
       y = "t-SNE Dimension 2") +
  theme_classic()
print(tsne_plot)
