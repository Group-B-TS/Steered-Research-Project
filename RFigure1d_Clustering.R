library(Rtsne)
library(ggplot2)

data <- read.csv("GSE202695_counts_afterQC.csv")
metadata <- read.csv("GSE202695_metadata.csv")

combined_data <- cbind(data, Model = metadata$model)

numeric_data <- combined_data[, sapply(combined_data, is.numeric)]

tsne_result <- Rtsne(numeric_data, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500, partial_pca=TRUE)

tsne_df <- data.frame(X1 = tsne_result$Y[,1], X2 = tsne_result$Y[,2],Model = combined_data$Model)

ggplot(tsne_df, aes(x =X1, y = X2, color = Model)) + geom_point() + labs(title = "tSNE Plot of Sequenced Cells by Model of Origin")