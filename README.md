
# Group B-TS Steered Research Project




## Table of Content

    Introduction
    Usage
    Authors
    Contributing


# Introduction 
This project aims to replicate the computational analysis outlined in the Hamelin et al. Single-cell Analysis Reveals Inter- and Intratumour Heterogeneity in Metastatic Breast Cancer paper, providing a comprehensive suite of scripts necessary to recreate their workflow. For researchers and computational biologists interested in exploring and reproducing the findings of the aforementioned paper, this repository contains a collection of scripts designed to facilitate various tasks such as differential gene expression analysis, clustering, quality control, and data visualization.


## Usage

This repository provides a set of scripts essential for both reanalysis of existing data and executing a novel pipeline. Below is a summary of each script and its function within the respective pipelines:

<details>
<summary>Reanalysis Pipeline</summary>

- Retrieve_SRA_Accessions.sh: Retrieves SRA accessions for dataset acquisition.
- download_script.sh: Downloads data from the retrieved SRA accessions.
- checksum_integrity.sh, checksum_verification.sh: Performs checksum integrity checks on downloaded data.
- fastqc.sh: Executes FastQC for quality assessment of raw sequencing data.
- multiqc.sh: Generates a consolidated report using MultiQC, summarizing results from various quality control analyses.
- fastqscreen.sh: Performs quality checks and screens raw sequencing data.
- STAR_alignment.sh: Aligns the raw sequencing reads to the reference genome using STAR.
- STAR_GenomeIndexing.sh: Indexes the reference genome for subsequent alignment steps.
- sort_index_samtools.sh: Sorts and indexes alignment files using Samtools.
- raw_counts.R: Generates raw count matrices from aligned reads.
- fastqscreen_aggregate.R: Aggregates FastQScreen results from both strands.
- filter_raw_counts.py: Filters raw count matrices.
- countsafterQC.csv: Resulting count matrices after quality control.
- Initial Version of Clustering Analysis: Implements the initial version of clustering analysis.
- CellCycleCode.R: Performs cell cycle analysis.
- DGEA.R: Executes differential gene expression analysis and gene set enrichment.
</details>

<details>
<summary>Novel Pipeline</summary>

- fastp_adapter_trim_pipeline2.sh: Trims adapter sequences from raw sequencing data using Fastp.
- hisat_alignment.sh: Aligns the trimmed reads to the reference genome using HISAT.
- build_hgfm_transcript_index.sh: Builds a transcript index for HISAT alignment.
- htseq_count.sh: Generates count matrices using HTSeq.
- Initial Clustering Novel Pipeline: Implements the initial clustering step for the novel pipeline.
- cell_cycle_clustering_gsea_novel_pipeline.R: Performs cell cycle correction and gene set enrichment analysis for the novel pipeline.
  
</details>


## Authors

-[@Nayah18](https://www.github.com/Nayah18)
-[@eleftheriosserafim](https://www.github.com/eleftheriosserafim)
-[@sheron-mathew](https://www.github.com/sheron-mathew)
-[@Alerocodes](https://www.github.com/Alerocodes)
-[@zaikanu99](https://www.github.com/zaikanu99)
-[@ikraanm](https://www.github.com/ikraanm)

## Contributing

Contributions are always welcome! If you're interested in fixing a bug, adding new features, or improving documentation, your contributions help make this project better for everyone.

