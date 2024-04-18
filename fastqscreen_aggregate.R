library(dplyr)

# Helper function to count the lines in a file
count_lines <- function(file) {
  length(readLines(file))
}

# Function to extract uniquely mapped counts for Human and Mouse from each file
extract_counts <- function(file) {
  # Read the file, assuming that the first row is the header, and skipping the footer line
  data <- read.table(file, header = TRUE, sep = "\t", skip = 1, fill = TRUE, strip.white = TRUE, comment.char = "", nrows = count_lines(file) - 2)
  human_counts <- sum(data[data$Genome == "Human", "X.One_hit_one_genome"], na.rm = TRUE)
  mouse_counts <- sum(data[data$Genome == "Mouse", "X.One_hit_one_genome"], na.rm = TRUE)
  id <- gsub("(_[12])?_screen\\.txt$", "", basename(file))  # Extract ID without the suffix
  
  return(data.frame(ID = id, Human = human_counts, Mouse = mouse_counts))
}

# List all files in the specified directory that match the pattern
files <- list.files(path = "/home/ajo26/SRP/fastqscreen_test", pattern = "screen\\.txt$", full.names = TRUE)

# Process each file to extract counts and bind rows together
results_df <- do.call(rbind, lapply(files, extract_counts))

# Aggregate the results by unique SRR run ID
final_results <- aggregate(cbind(Human, Mouse) ~ ID, data = results_df, FUN = sum)

# Output the results to a CSV file
write.csv(final_results, "fastqscreen_counts.csv", row.names = FALSE)

print(final_results)
