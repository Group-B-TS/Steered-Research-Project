#!/usr/bin/bash
#SBATCH --job-name=fastq_screen test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ajo26@student.le.ac.uk
#SBATCH --export=NONE

# Full path to the FastQ Screen executable
FASTQ_SCREEN_EXEC=/home/a/ajo26/fastq_screen/FastQ-Screen-0.15.3/fastq_screen

# Directory where your FASTQ files are located
FASTQ_DIR=/home/a/ajo26/fastq_screen/FastQ-Screen-0.15.3/fqs_test_dataset.fastq.gz 

# Directory to output the FastQ Screen results
OUTPUT_DIR=/home/a/ajo26/fastq_screen/FastQ-Screen-0.15.3/fastqscreen_results


# FastQ Screen configuration file
CONFIG_FILE=/home/a/ajo26/fastq_screen/FastQ-Screen-0.15.3/fastq_screen.conf

# Loop through all FASTQ files in the directory
for fastq_file in ${FASTQ_DIR}/*.fastq
do
    echo "Processing file ${fastq_file}"
    # Run FastQ Screen for each FASTQ file
    ${FASTQ_SCREEN_EXEC} --conf ${CONFIG_FILE} --outdir ${OUTPUT_DIR} ${fastq_file}
done

echo "FastQ Screen processing completed."
