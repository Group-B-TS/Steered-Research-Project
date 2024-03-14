## Steered-Research-Project shell script for download of raw data and conversion into fastq file ##

#!/bin/bash

TSV_FILE="/scratch/alice/i/ib156/Steered_project/PRJNA837007-info.tsv"

while IFS=$'\t' read -r srx srr platform
do
    echo "SRR ID: $srr"  #checking which srr ID is being downloaded, was used for debugging
    echo "Downloading $srr.."
    prefetch $srr                #download using srr accession and using sra toolkit prefetch command
    fastq-dump --split-files $srr #conversion into fastq, using sra toolkit fastq-dump command
done < "$TSV_FILE"

echo "All files were downloaded."
