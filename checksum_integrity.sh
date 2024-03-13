#!/bin/bash

#SBATCH --job-name=checksum_fastq
#SBATCH --output=checksum_fastq_%j.out
#SBATCH --error=checksum_fastq_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --partition=medium
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk

# Checking if the fastq_checksums.md5 file already exists and removing it to avoid appending to an old file (due to trial/errors)
if [ -f fastq_checksums.md5 ]; then
    rm fastq_checksums.md5
fi

# Generating checksums for all FASTQ files starting with SRR and ending with .fastq
for file in SRR*.fastq; do
    md5sum "$file" >> fastq_checksums.md5
done
