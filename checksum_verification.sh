#!/bin/bash

#SBATCH --job-name=checksum_verification
#SBATCH --output=checksum_verification_%j.out
#SBATCH --error=checksum_verification_%j.err
#SBATCH --time=3-00:00:00 #much more time than needed just in case due to number of files
#SBATCH --ntasks=1
#SBATCH --mem=8G #likely more than enough as well
#SBATCH --partition=medium
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk

cd /scratch/alice/i/ib156/Steered_project/
md5sum -c fastq_checksums.md5 #due to the HPC system being temporarily unavailable, we executed this checksum verification command directly within our scratch directory
