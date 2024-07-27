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
### 11. cannaspy-main-project
### 12. communicating_with_synapse-master-project
### 13. datarail-master-project
### 14. mcmicro-master-project
### 15. scRNASeq-pilot project
### 16. single_cell_RNAseq_Visualization-main-project
### 17. workflows-master-project
