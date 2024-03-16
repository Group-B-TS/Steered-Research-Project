#!/bin/bash

#SBATCH --job-name=multiqc_summary
#SBATCH --output=multiqc_%j.out
#SBATCH --error=multiqc_%j.err
#SBATCH --time=1-00:00:00  #More than enough time to complete the task...
#SBATCH --mem=4G
#SBATCH --partition=medium
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk

#Multiqc module version on Alice
module load py-multiqc/1.14-ateq76h

# Navigating to the directory containing our FastQC results (just in case)
cd /scratch/alice/i/ib156/Steered_project/fastqc_results

# Running MultiQC to aggregate the FastQC reports into a summary
multiqc . -o .
