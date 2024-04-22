# Loading necessary libraries
library(GenomicRanges)
library(ensembldb)
library(bamsignals)
library(Rsamtools)  # Rsamtools included for the handling of BAM files

# Loading Ensembl database
ensDbPath <- "/scratch/alice/i/ib156/Steered_project/EnsDb_Homo_sapiens_GRCh38_ensembl_96.sqlite"
ensDb <- EnsDb(ensDbPath)

# Geting chromosome lengths and preparing them (to have same naming convention across files, as recommended by bamcount guide)
chr_lengths <- seqlengths(ensDb)
names(chr_lengths) <- paste0("chr", names(chr_lengths))

desired_seqlevels <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                       "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
                       "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
chr_lengths <- chr_lengths[names(chr_lengths) %in% desired_seqlevels]

# Generating tiles of 20 bp across the genome
tiles <- tileGenome(chr_lengths, tilewidth = 20, cut.last.tile.in.chrom = TRUE)

# Confirming that tiles is a GenomicRanges object, was useful for debugging
if (!inherits(tiles, "GenomicRanges")) {
  stop("Tiles object is not a GenomicRanges object.")
}

# Defining the directory containing BAM files and the pattern for bam files
bam_dir <- "/scratch/alice/i/ib156/Steered_project/sorted_bam_files"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)


# Directory to save individual coverage CSV files
coverage_output_dir <- "/scratch/alice/i/ib156/Steered_project/coverage_files/"

# Processing each BAM file, was useful for debugging
for (bam_file in bam_files) {
  cat("Processing file:", bam_file, "\n")
  
  tryCatch({    # Useful for debugging too
    coverage <- bamCount(bam_file, gr = tiles, mapqual = 0, shift = 0, ss = TRUE, paired.end = "filter", tlenFilter = c(0, 1000), filteredFlag = -1, verbose = TRUE)
    if (!is.null(coverage)) {
      # Generating a unique filename for the coverage data
      coverage_filename <- paste0(coverage_output_dir, basename(bam_file), "_coverage.csv")
      write.csv(coverage, coverage_filename, row.names = FALSE)
      cat("Coverage calculated and saved successfully for file: ", bam_file, "\n")
    } # Initially, crash and nonsensical errors occured after looping through multiple files, then after a single file...
  }, error = function(e) {
    cat("Error processing file:", bam_file, " - Error: ", e$message, "\n")  # No error message is displayed before R crashes
  })
}

# After processing all files, we combine all individual coverage files
coverage_files <- list.files(coverage_output_dir, pattern = "_coverage.csv$", full.names = TRUE)
combined_coverage <- do.call("rbind", lapply(coverage_files, read.csv, header = TRUE))

# Writing the combined data frame to a new CSV file, like the authors
final_coverage_file <- paste0(coverage_output_dir, "final_count_matrix.csv")
write.csv(combined_coverage, final_coverage_file, row.names = TRUE)
cat("All coverage files have been merged into:", final_coverage_file, "\n") 

# Loading GTF file and convert to SAF format (necessary according errors)
gtf_path <- "/home/ib156/Desktop/Homo_sapiens.GRCh38.96_chr.gtf"
gtf_data <- import(gtf_path)
saf_data <- data.frame(
  GeneID = mcols(gtf_data)$gene_id,
  Chr = seqnames(gtf_data),
  Start = start(gtf_data),
  End = end(gtf_data),
  Strand = ifelse(strand(gtf_data) == "+", 1, "-1")
)

# Definining SAF file path
saf_path <- "/home/ib156/Desktop/Homo_sapiens.GRCh38.96.saf"
write.table(saf_data, saf_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Running featureCounts on BAM files and saving individual count CSV files
for (bam_file in bam_files) {
  output_file <- paste0(counts_output_dir, basename(bam_file), "_counts.csv")
  counts <- featureCounts(files = bam_file,
                          annot.ext = saf_path,
                          isPairedEnd = TRUE,
                          strandSpecific = 1,
                          countChimericFragments = TRUE)
  counts_df <- as.data.frame(counts$counts, check.rows = FALSE, check.names = FALSE)
  rownames(counts_df) <- counts$annotation$GeneID
  write.csv(counts_df, output_file, row.names = TRUE)
  cat("Completed featureCounts for", bam_file, "\n")
}

# Function to read CSV files with gene IDs as row names
read_counts <- function(file) {
  read.csv(file, row.names = 1, check.names = FALSE)
}

# Reading and combining count CSV files
count_files <- list.files(counts_output_dir, pattern = "_counts\\.csv$", full.names = TRUE)
all_counts_list <- lapply(count_files, read_counts)

# Verifying consistency in row names
if (!all(sapply(all_counts_list, function(df) all(rownames(df) == rownames(all_counts_list[[1]]))))) {
  stop("Mismatch in row names/order among count files.")
}

# Combining count data frames with cbind
all_counts <- do.call(cbind, all_counts_list)

# Writing the combined counts to the final CSV file
final_counts_file <- paste0(counts_output_dir, "final_count_matrix.csv")
write.csv(all_counts, final_counts_file, row.names = TRUE)
cat("All count files have been merged into:", final_counts_file, "\n")
