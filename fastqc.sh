#!/bin/bash

#SBATCH --job-name=fastqc_analysis
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=medium
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk

# Loading FastQC module
module load fastqc/0.12.1-hkgpcde

# Defining the directory for FastQC results
RESULT_DIR="/scratch/alice/i/ib156/Steered_project/fastqc_results"

# Ensuring the result directory exists, if not, it will be created
mkdir -p "$RESULT_DIR"


# Changing directory  in case it is not run from the right one
cd /scratch/alice/i/ib156/Steered_project/

# Running FastQC on all FASTQ files starting with SRR and ending with .fastq
# Saving the results to the above-mentioned results directory
fastqc -t 4 SRR*.fastq -o "$RESULT_DIR"
