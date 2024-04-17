library(Rtsne)
library(ggplot2)

setwd("/alice-home/1/e/es452/Research_Project/")

# Read in the data
data <- read.csv("GSE202695_counts_afterQC.csv", header = TRUE)
data <- data[,-1]  # Assuming the first column is not part of the data

meta <- read.csv("GSE202695_metadata.csv", header = TRUE)

# Log-normalize the data by adding 1 to avoid log(0) issues
data_log_normalized <- log2(data + 1)

# Transpose the data so that samples are rows and genes are columns for PCA
data_t_log_normalized <- t(data_log_normalized)

# Identify and remove constant/zero-variance columns before PCA
non_constant_columns <- apply(data_t_log_normalized, 2, var) != 0
data_t_log_normalized <- data_t_log_normalized[, non_constant_columns]

pca_results <- prcomp(data_t_log_normalized, scale. = FALSE)

# Choose a number of principal components based on variance explained
pca_data <- pca_results$x[, 1:50]

# Set the seed for reproducibility before running t-SNE
set.seed(42)
tsne_results <- Rtsne(pca_data, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000)

# Convert t-SNE results to a data frame
tsne_data <- as.data.frame(tsne_results$Y)

# Ensure that the rownames from your PCA match your metadata
rownames(meta) <- rownames(pca_data)

# Combine the t-SNE results with the metadata
tsne_data <- cbind(tsne_data, meta)

# Define the colours for the new names
name_mapping <- c("HBRX1921" = "PDX1",
                  "HBRX2344" = "PDX4",
                  "HBRX2353" = "PDX2",
                  "HBRX3078" = "PDX3",
                  "MDAMB231" = "MDAMB231")

# Map the 'model' column in meta to the new names using the mapping vector
meta$Model <- name_mapping[meta$model]

# Define the colors for the new names
my_colors <- c("PDX3" = "red", "PDX4" = "blue", "PDX1" = "green", "MDAMB231" = "orange", "PDX2" = "purple")

tsne_data$Model <- meta$Model

tsne_plot <- ggplot(tsne_data, aes(x = V1, y = V2, color = Model)) +
  geom_point() +
  labs(title = "t-SNE Plot of Sequenced Cells by Model of Origin",
       x = "TSNE.logcounts 1",
       y = "TSNE.logcounts 2") +
  theme_classic() 
print(tsne_plot)
