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
