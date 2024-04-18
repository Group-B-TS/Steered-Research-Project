#!/bin/bash

#SBATCH --job-name=fastp_trimming
#SBATCH --output=fastp_output_%j.txt
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=8G


# Activating the Conda environment
source /home/i/ib156/miniconda3/etc/profile.d/conda.sh
conda activate fastp_env


# Defining our directory paths
input_dir="/scratch/alice/i/ib156/Steered_project/"
output_dir="/scratch/alice/i/ib156/Steered_project/trimmed_fastq/"

# Ensuring output directory exists
mkdir -p "$output_dir"

# Looping through all files that match SRR*_1.fastq in the input directory
for file1 in "${input_dir}"SRR*_1.fastq; do
    file2="${file1/_1.fastq/_2.fastq}"  # Constructing the name of the mate pair file

    # Outputting filenames based on input filenames
    out1="${output_dir}/$(basename "${file1/_1.fastq/_1.trimmed.fastq}")"
    out2="${output_dir}/$(basename "${file2/_2.fastq/_2.trimmed.fastq}")"
    html_report="${output_dir}/$(basename "${file1/_1.fastq/.report.html}")"
    json_report="${output_dir}/$(basename "${file1/_1.fastq/.report.json}")"

    # Running fastp
    fastp -i "$file1" -I "$file2" -o "$out1" -O "$out2" --html "$html_report" --json "$json_report"
done

# Deactivating our Conda environment
conda deactivate

echo "Fastp processing complete." # Not necessary for SLURM submission, but useful for local runs


