#!/bin/bash

# Defining  the directory containing our BAM and BAI files
BAM_DIR="/scratch/alice/s/sm1073/steered_project2/"

# Defining the directory to save the log files
LOG_DIR="/scratch/alice/i/ib156/Steered_project/logs"
mkdir -p "${LOG_DIR}"

# Creating the two log files for corrupted files
CORRUPTED_BAM_LOG="${LOG_DIR}/corrupted_bam_files.log"
CORRUPTED_BAI_LOG="${LOG_DIR}/corrupted_bai_files.log"

# Clearing the log files just in case
> "${CORRUPTED_BAM_LOG}"
> "${CORRUPTED_BAI_LOG}"

# Checking the BAM and BAI files
for bam_file in "${BAM_DIR}"/*.bam; do
  # Checking BAM file integrity
  if ! samtools quickcheck "$bam_file"; then
    echo "Corrupted BAM: $bam_file" >> "${CORRUPTED_BAM_LOG}"
  fi

  # Constructing the BAI filename
  bai_file="${bam_file}.bai"

  # Checking if the BAI file exists
  if [ ! -f "$bai_file" ]; then
    echo "Missing BAI: $bai_file" >> "${CORRUPTED_BAI_LOG}"
    continue
  fi

  # Checking BAM index (BAI) integrity by attempting to use it
  if ! samtools idxstats "$bam_file" > /dev/null; then
    echo "Corrupted or Incompatible BAI: $bai_file" >> "${CORRUPTED_BAI_LOG}"
  fi
done

# Outputting the result
echo "BAM file integrity check complete. Corrupted files listed in ${CORRUPTED_BAM_LOG}"
echo "BAI file integrity check complete. Corrupted or missing files listed in ${CORRUPTED_BAI_LOG}"
