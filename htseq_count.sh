#!/bin/bash

#SBATCH --job-name=htseq_count   
#SBATCH --cpus-per-task=32                 
#SBATCH --mem=32G                        
#SBATCH --time=72:00:00                  
#SBATCH --output=htseq_count_%j.log       
#SBATCH --error=htseq_count_%j.err        
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ajo26@student.le.ac.uk

# Path to the directory containing sorted BAM files
bam_dir="/scratch/alice/a/ajo26/SRP/sorted_bam"
# Path to the GTF file
gtf_file="/scratch/alice/a/ajo26/SRP/hisat2/genome.gtf"

# Path to the output directory where you want to store count files
output_dir="/scratch/alice/a/ajo26/SRP/count_dir"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Navigate to the directory with BAM files
cd "$bam_dir"

# Loop through each sorted BAM file in the directory with the specific naming pattern
for bamfile in *sorted_output.bam; do
    # Define output filename by replacing 'sorted_output.bam' with 'counts.txt'
    # and prefixing with the output directory path
    output_file="${output_dir}/${bamfile/%sorted_output.bam/counts.txt}"
    
    echo "Processing $bamfile"
    
    # Run htseq-count with the specified options
    htseq-count --mode=union \
                --stranded=yes \
                --format=bam \
                --order=pos \
                --idattr=gene_id \
                -a 10 \
                -t exon \
                "$bamfile" \
                "$gtf_file" > "$output_file"
    
    echo "Output written to $output_file"
done

