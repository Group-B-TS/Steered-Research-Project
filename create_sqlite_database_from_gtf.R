library(ensembldb)

# Defining file paths
gtfPath <- "/home/ib156/Desktop/Homo_sapiens.GRCh38.96.gtf"
ensDbPath <- "/home/ib156/Desktop/EnsDb_Homo_sapiens_GRCh38_ensembl_96.sqlite"

# Creating the EnsDb database from the GTF file
ensDbFromGtf(gtf = gtfPath, 
             outfile = ensDbPath, 
             organism = "Homo sapiens", 
             genomeVersion = "GRCh38", 
             version = "96")

# Loading the EnsDb object directly
edb <- EnsDb(ensDbPath)

# Fetching sequence information including chromosome lengths
seqinfo <- seqinfo(edb)

# Printing the seqinfo object to see details of each chromosome/contig
print(seqinfo)

# Calculating and printing the total genome size for potential future use...
genomeSize <- sum(seqlengths(seqinfo))
cat("Total genome size:", genomeSize, "bases\n")

