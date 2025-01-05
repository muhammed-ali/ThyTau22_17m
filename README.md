# Longitudinal Single-Cell Transcriptomic Analysis of Sex-Dependent Changes from THY-Tau22 mouse model of AD

# Table of contents
* [Introduction](#introduction)
* [Content](#content)
* [Data](#data)
* [Requirements](#requirements)
* [License](#license)
* [Instructions](#instructions)

# Introduction
This repository contains the code for bioinformatics analyses described in the article "Longitudinal single-cell transcriptomic analysis reveals sex-dependent changes in THY-Tau22 mice at early and late stages of tau pathology".

This project investigated sex-dependent molecular changes in the THY-Tau22 mouse model of Alzheimer's disease through single-cell RNA sequencing analysis, comparing transcriptomic profiles at 7 and 17 months of age to understand disease progression and sex-specific responses.

# Content
The code covers the following main analysis steps:

1. QC, Normalization (SCT), Clustering, and Cell Type Annotation
2. Differential Expression Analysis (bulk and cell type-specific)
3. DEA using edgeR for sex-interaction analysis  
4. Pathway Enrichment Analysis
5. Gene Regulatory Network (GRN) Analysis 
6. Analysis of DEG Overlaps between Datasets
7. Longitudinal Analysis (7 to 17 months temporal changes)
   
# Data
The single-cell RNA sequencing data is available in the NCBI Gene Expression Omnibus (GEO) database:
- THY-Tau22 mice (7 months): GSE245035
- THY-Tau22 mice (17 months): GSE285506
- Tg2576 mouse model: [accession pending]
- Human AD cortical tissue: GSE138852

# Requirements
The code was written in R (version 4.2.2) and relies on multiple R and Bioconductor packages, including:
- Seurat (v4.3.0)
- clusterProfiler 
- enrichplot
- cowplot
- ggplot2
- HGNChelper
- ggVennDiagram
- Additional packages listed at the beginning of each R script

# License
The code is available under the MIT License.

# Instructions
The code was tested on R 4.2.2 on both current Mac and Linux operating systems, but should be compatible with later versions of R installed on current Mac, Linux or Windows systems.

Required R packages can be installed using:

For CRAN packages:

    install.packages("BiocManager::install("ADD_NAME_OF_THE_PACKAGE")

R packages from Bioconductor can be installed with the following commands:

    if (!require("BiocManager", quietly = TRUE))

        install.packages("BiocManager")

    BiocManager::install("ADD_NAME_OF_THE_PACKAGE")

To run the code, the correct working directory containing the input data must be specified at the beginning of the R-scripts, otherwise the scripts can be run as-is.

The scripts should be run in the following order:

    quality_control_clustering.R

    differential_expression.R

    sex_interaction_analysis.R

    pathway_enrichment.R

    grn_analysis.R

    deg_overlap_analysis.R

    longitudinal_analysis.R
