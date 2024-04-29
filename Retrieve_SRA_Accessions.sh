#Commands to get all SRA accessions related to the author's PRJNA accession in a tsv file 
#Adapted from: https://hackmd.io/@AstrobioMike/SRA-fastq-download#Downloading-fastq-files-based-on-a-bioproject-PRJNA…
#Requires prior installation of the NCBI Entrez Direct (EDirect) utilities

esearch -db sra -query PRJNA837007 | esummary | xtract -pattern DocumentSummary -element Experiment@acc,Run@acc,Platform@instrument_model > PRJNA837007-info.tsv

head PRJNA837007-info.tsv
