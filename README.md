# Immunitas
Github Examples for Sr. Computational Biologist position 

# Project Organization

I have been involved with several organizations and projects, each with its own repository. To streamline this, I have consolidated these projects into a single GitHub repository. The selected examples span various fields, including single-cell research, immunology, oncology, and data processing. To maintain confidentiality, I have anonymized company names except for H3 Biomedicine, which has since ceased operations. This repository highlights 17 projects in total, covering a broad scope of my work. While some projects were solely my responsibility as a coder or analyst, others were collaborative efforts, resulting in variations in format and style.

# Project Organization

### 1. BCMA-SMLA-project

Cancer cells commonly develop resistance to immunotherapy by loss of antigen expression. Here, we use our CRISPR interference– and CRISPR activation–based functional genomics platform to systematically identify pathways controlling cell surface expression of the multiple myeloma immunotherapy antigen B-cell maturation antigen (BCMA). We discovered that pharmacologic inhibition of HDAC7 and the Sec61 complex increased cell surface BCMA. This contains the analysis that led to discovering those pathways and targets that led to an expanded project. 

### 2. CyCif_Manager-master-project

The CyCIF Manager project aims to provide a streamlined pipeline platform for CyCIF (Cyclic Immunofluorescence) analysis, both on local machines and on the O2 cluster at Harvard Medical School (HMS). This project includes infrastructure for managing and analyzing high-dimensional imaging data, essential for research in fields like immunology and oncology. The platform ensures consistent and efficient processing of CyCIF data, facilitating the integration of different datasets and enhancing reproducibility in analysis.

Key features of the CyCIF Manager include scripts and tools for setting up the necessary environments, managing pipeline versions, and processing example datasets. The repository also includes detailed documentation and installation guides to help users deploy and utilize the pipeline effectively.

### 3. DrugResponse-master-project

The DrugResponse project, hosted on GitHub, is an R package designed for predicting the response of cancer cell lines or patient samples to seven FDA-approved cancer drugs: Gemcitabine, Gefitinib, Cisplatin, Doxorubicin, Docetaxel, Paclitaxel, and Carboplatin. This package uses gene expression data from Affymetrix U133 Plus 2.0 arrays to make these predictions.

Key Features:
Prediction Function: The main function of the package, DrugResponse.predict, takes gene expression data and predicts the likelihood of a positive or negative response to each of the seven drugs.
Data Requirements: The input data can be in the form of CEL files (raw microarray data) or preprocessed gene expression files.
Pipeline Scripts: The project includes bash scripts for running tests on sample data, making it easier to process multiple samples in batch mode.

### 4. EPAT-project
This document is for analyzing EPAT - Chempartner data. The markdown file would do correlation of Relative IC50 values to mRNA and mutation data from CCLE cell-lines data.

The EPAT (Electronic Phenotype and Annotation Tool) project on GitHub is a comprehensive tool designed for annotating genetic variants. The project provides scripts for analyzing VCF (Variant Call Format) files against reference genomes and annotations. The core functionality includes assessing the impact of genetic mutations on protein function, using tools such as PROVEAN (Protein Variation Effect Analyzer) to predict whether mutations are harmful.

Key features of the EPAT project include:

Input Requirements: Users need to provide a VCF file for input, a FASTA file for the reference genome, and a GTF file for annotation.
Docker and Singularity Support: The tool supports running analyses through Docker and Singularity containers, ensuring compatibility and ease of setup across different environments.
Output: The analysis results in a detailed output file that includes columns for PROVEAN scores and predictions, indicating the potential impact of each mutation on protein function.

### 5. SF3B1_WestLake_CN-project

The "SF3B1_WestLake_CN-project" is focused on analyzing the impact of mutations in the SF3B1 gene, particularly in relation to splicing factor 3b, subunit 1 (SF3B1) mutations. This project involves extensive genomic and bioinformatic analyses to understand the mutations' roles in various cancers, specifically chronic lymphocytic leukemia (CLL). The project utilizes various bioinformatics tools and methods to analyze sequencing data, generate insights into mutation effects, and potentially identify therapeutic targets.

### 6. SRA_Processor-main-project

This data pipeline is focused on the ability to automate a parallized searching, downloading, and formatting sequencing data based on a string of search words for future pipeline work.

### 7. TAM-PARP-2019-master-project

The TAM-PARP-2019 project focuses on overcoming PARP-inhibitor resistance in breast cancer by harnessing tumor-associated macrophages. The repository contains code for processing and analyzing CyCif and RNA-Seq data, with scripts for data stitching, segmentation, and gene set enrichment analysis. Key technologies used include Ashlar for image stitching, Sleuth for differential RNA-seq analysis, and fGSEA for gene set enrichment, supporting the research detailed in their associated publication.

### 8. [company name]-nanostring-analyses

This project contains several projects that were part of a much larger effort to manage data analysis efforts for a company that develops patient-derived cancer organoids that pharmaceutical companies would contract to trial their drugs to test efficacy, mechanism of action studies, and general drug/target evaluation. The various evaluations largely utilized nanostring as well as image analysis to form the basis for their organoid platform such as pooling samples, sensitivity studies, as well as some case studies for known and unknown drug mechanism of action studies.

### 9. [company name]-scRNAseq-Analysis-main-project

This project is dedicated to single-cell RNA sequencing analysis that contains detailed scRNAseq analysis, focusing on advanced computational workflows to understand cellular dynamics and gene expression profiles at the single-cell level. The client sought to understand how their drug targets influenced various aspects on the mouse models they deployed.

### 10. [company name]_mc38_scRNAseq-main-project

This project was dedicated to a single-cell RNA sequencing analysis for a client that was evaluating several publically available datasets to gain insight into how they wish to proceed for their drug development with regards to targets.

### 11. cannaspy-main-project

This project was a data gathering project for a couple of websites to webscrape their data into a dedicated format for future data analysis purposes.

### 12. communicating_with_synapse-master-project

This project was designed to assist with automating of data formatting, syncing with a cloud based platform, and provide a general version controlled resource for data tracking.

### 13. datarail-master-project

This system is specifically designed to support the design, analysis, and visualization of high-throughput dose-response experiments. The project involves the use of various tools and frameworks to facilitate these experiments, aiming to improve the efficiency and effectiveness of data handling and analysis processes in experimental settings.

### 14. mcmicro-master-project

The "mcmicro-master-project" on GitHub, forked from the labsyspharm/mcmicro repository, is an extensive multiple-choice microscopy pipeline developed for processing multiplexed whole slide imaging and tissue microarrays. This end-to-end pipeline includes steps like stitching and registration, segmentation, and single-cell feature extraction. Each phase of the pipeline is designed to be containerized, enhancing portability across different computing environments.

This comprehensive system allows researchers to handle complex imdaging data systematically, facilitating detailed analysis at the cellular level, which is critical in areas such as pathology and cancer research. The pipeline is supported by various grants and is part of collaborative efforts that underline its significance in advancing medical and biological research.

### 15. scRNASeq-pilot project

This project is for a client that was working on establishing best practices for their single cell RNA-Seq projects for drug discovery applications. This analysis compared a couple of different vendors as well as as the pro's and con's on different analysis pipelines to gain understanding into their various drug discovery projects.

### 16. single_cell_RNAseq_Visualization-main-project

This project was for a client that sought different single cell RNA-Seq visualization perspectives to gain insight into their drug discovery projects.

### 17. workflows-master-project

This project purpose was to streamline and automate various workflows in order how to structure and execute bioinformatics workflows efficiently, ensuring scalability and reproducibility. This project used Nextflow, a workflow tool that simplifies the deployment and sharing of complex pipeline projects across different computational environments. This project incorporates a series of custom scripts (BASH, R, Python, etc.) that are managed consistently within the Nextflow framework, demonstrating best practices for building and maintaining reproducible scientific workflows. Key features of this project include 1) Use of Docker containers to ensure consistency across different computational environments, 2) Integration with batch schedulers like SLURM for efficient job management on computational clusters, 3) Configuration profiles allowing different settings for different execution environments, such as local or cloud-based systems.

