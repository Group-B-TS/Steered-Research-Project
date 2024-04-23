#!/bin/bash

#SBATCH --job-name=hisat_alignment   
#SBATCH --cpus-per-task=32                 
#SBATCH --mem=32G                        
#SBATCH --time=48:00:00                   
#SBATCH --output=hisat_alignment_%j.log       
#SBATCH --error=hisat_alignment_%j.err        
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ajo26@student.le.ac.uk

# Load the SAMtools module
module load samtools/1.18

# Activate Conda environment with HISAT2
source ~/miniconda3/bin/activate 
conda activate base

# Define paths
FASTQ_DIR="/scratch/alice/i/ib156/Steered_project/trimmed_fastq"
GENOME_INDEX="/scratch/alice/a/ajo26/SRP/hisat2"
OUTPUT_DIR="/scratch/alice/a/ajo26/SRP/sorted_bam" 

# Creating output directory 
mkdir -p $OUTPUT_DIR

# Navigate to the directory with FASTQ files
cd $FASTQ_DIR

# Loop through all paired-end read files
for file in *_1.trimmed.fastq; do
    sample_name=$(basename $file _1.trimmed.fastq)
    echo "Processing $sample_name"

    # Alignment, conversion to BAM, sorting, and indexing in one go
    hisat2 -p 32 -x $GENOME_INDEX/genome_tran -1 ${file} -2 ${sample_name}_2.trimmed.fastq | \
    samtools view -bS - | \
    samtools sort -@ 32 -m 4G -o $OUTPUT_DIR/${sample_name}_sorted_output.bam - && \
    samtools index $OUTPUT
done
