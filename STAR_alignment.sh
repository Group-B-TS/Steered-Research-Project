#!/bin/bash

#SBATCH --job-name=star_alignment
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4-00:00:00 #plenty of time
#SBATCH --partition=medium
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk

module load star/2.7.10b-m3zkpic

GENOME_DIR="/scratch/alice/i/ib156/Steered_project/indexing"
FASTQ_DIR="/scratch/alice/i/ib156/Steered_project"
OUTPUT_DIR="/scratch/alice/i/ib156/Steered_project/alignments"
mkdir -p "$OUTPUT_DIR"    # Ensuring the output directory exists

for R1 in ${FASTQ_DIR}/*_1.fastq; do    # Looping through all FASTQ files ending with '_1.fastq' to identify forward reads
    R2=${R1%_1.fastq}_2.fastq      # Generating the corresponding reverse read file name by replacing '_1.fastq' with '_2.fastq'
    SAMPLE=$(basename ${R1} _1.fastq)  # Extracting the sample name by removing the '_1.fastq' suffix
    # Running STAR aligner for each pair of forward and reverse reads
    STAR --genomeDir $GENOME_DIR \
         --readFilesIn $R1 $R2 \
         --runThreadN $SLURM_CPUS_PER_TASK \
         --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE}_  #Output files prefixes are based on the corresponding sample name in the specified output directory
done
