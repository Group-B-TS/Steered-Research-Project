setwd("/alice-home/1/e/es452/Research_Project/")

library(Rtsne)
library(ggplot2)

data <- read.csv("GSE202695_counts_afterQC.csv", header = TRUE)
data <-data[,-1]

meta <- read.csv("GSE202695_metadata.csv", header = TRUE)

# Transpose data so that samples are rows and genes are columns
data_t <- t(data)

# Convert to data frame
data_t <- as.data.frame(data_t)

# Match the rows of data_t with rows of meta
# Assume first column of meta is the identifier
rownames(data_t) <- meta[, 1]

# Perform PCA for dimensionality reduction
pca_results <- prcomp(t(data_t), scale. = TRUE)  

# Choose a number of principal components to keep
pca_data <- as.data.frame(pca_results$x[,1:50])

# Run t-SNE on the PCA-reduced data
set.seed(40)
tsne_results <- Rtsne(as.matrix(data_t), dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

# Combine the t-SNE output with the categories from metadata
tsne_data <- as.data.frame(tsne_results$Y)
tsne_data$Model <- meta[, 'model']  

# Plotting the t-SNE results colored by the extracted categories
tsne_plot <- ggplot(tsne_data, aes(x = V1, y = V2, color = Model)) +
  geom_point() +
  labs(title = "t-SNE Plot of Sequenced Cells by Model of Origin",
       x = "TSNE.logcounts 1",
       y = "TSNE.logcounts 2") +
  theme_classic()
print(tsne_plot)
