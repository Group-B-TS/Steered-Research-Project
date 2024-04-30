#!/bin/bash

#SBATCH --job-name=hisat2_index_build    
#SBATCH --cpus-per-task=16                 
#SBATCH --mem=160G                         
#SBATCH --time=24:00:00                    
#SBATCH --output=hisat2_index_%j.log       
#SBATCH --error=hisat2_index_%j.err        
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ajo26@student.le.ac.uk
                

# Activate Conda environment with HISAT2
source ~/miniconda3/bin/activate 
conda activate base

# Set the directory where all the data and index files are stored
DATA_DIR=/scratch/alice/a/ajo26/SRP/hisat2

cd $DATA_DIR

# Assuming genome.fa and genome.gtf are already in $DATA_DIR

# Extract splice sites and exons from the GTF file
hisat2_extract_splice_sites.py genome.gtf > genome.ss
hisat2_extract_exons.py genome.gtf > genome.exon

# Build the HGFM index with transcripts
hisat2-build -p 16 --exon genome.exon --ss genome.ss genome.fa genome_tran

echo "Index build complete. Files generated:"
