#!/bin/bash

#SBATCH --job-name=SortIndexSamtools
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00  
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk

# Loading Samtools module available on Alice
module load samtools/1.18

# Defining the directory where SAM files are stored
DIR="/scratch/alice/i/ib156/Steered_project/alignments"

# Looping through all SAM files in the directory
for SAM_FILE in ${DIR}/*_Aligned.out.sam; do
  # Extracting the base of the file name for output naming
  BASE_NAME=$(basename ${SAM_FILE} _Aligned.out.sam)

  # Sorting and converting to BAM
  samtools sort -o ${DIR}/${BASE_NAME}_sorted.bam ${SAM_FILE}

  # Indexing the sorted BAM file
  samtools index ${DIR}/${BASE_NAME}_sorted.bam
done
