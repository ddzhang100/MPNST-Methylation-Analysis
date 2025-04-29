# MPNST/pNF Methylation Profiling

**Author:** Daniel Zhang  
**Date:** April 2025

## Overview

This repository contains R analysis scripts for DNA methylation profiling in Neurofibromatosis type 1 (NF1) associated benign plexiform neurofibromas (pNF) and malignant peripheral nerve sheath tumors (MPNST), using both the UCSF and Toronto cohorts.  We run parallel pipelines using:

- **Minfi** (Illumina EPIC arrays)  
- **Sesame** (IDAT‐level normalization)  

and compare clustering, PCA, heatmaps, and differential methylation.

## Prerequisites

- R ≥ 4.1  
- [Bioconductor](https://bioconductor.org/)  
- Git (for version control / pushing to GitHub)

## R Package Dependencies

Install via CRAN or Bioconductor:

```r
install.packages(c("ggplot2","pheatmap","cluster","factoextra"))
BiocManager::install(c(
  "minfi",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "sesame",
  "sesameData",
  "ConsensusClusterPlus"
))
