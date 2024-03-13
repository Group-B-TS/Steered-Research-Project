#SBATCH --job-name=genomeIndexing
#SBATCH --output=genomeIndexing_%j.out
#SBATCH --error=genomeIndexing_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=medium
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ib156@student.le.ac.uk

# Loading the STAR module available on alice
module load star/2.7.10b-m3zkpic

# Running STAR commands for genome indexing
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /scratch/alice/i/ib156/Steered_project/indexing \
     --genomeFastaFiles /scratch/alice/i/ib156/Steered_project/indexing/hg38.analysisSet.fa \
     --sjdbGTFfile /scratch/alice/i/ib156/Steered_project/indexing/hg38.refGene.gtf \
     --genomeSAindexNbases 12

