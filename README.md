# Multi-tissue Mendelian Randomization Analysis of Oral Cancer Susceptibility
This repository contains the complete analysis pipeline for the manuscript submitted to Research Connections.

## Requirements
R 4.x (see sessionInfo.txt for full package versions)
PLINK 1.9 (Windows 64‑bit)

## How to reproduce
Clone or download this repository.
Edit WORK_DIR in analysis_pipeline.R to your local path.
Download required raw data as described in data/README_data.txt.
Run source("analysis_pipeline.R") in RStudio.

After running the script, local folders `results/` and `results/figures/` will be created automatically to store output tables and figures. These output files are not uploaded to GitHub.

## Data availability
All raw data are publicly available (GTEx, FinnGen, GWAS Catalog, TCGA). Due to file size limits, raw data are not included here. See data/README_data.txt for download links.

## Results
Output tables and figures will be generated locally under `results/` and `results/figures/` after successful execution of the analysis script. Result files are not hosted in this repository.

## Contact
Corresponding author: lizhg@lzu.edu.cn

## MR Reporting Standard (STROBE-MR Checklist Compliance)

This two-sample multi-tissue MR analysis fully adheres to the STROBE-MR reporting guidelines. All analytical steps required by the checklist are explicitly implemented in `analysis_pipeline.R`, including:

1. SNP instrument selection (P-value threshold, LD clumping parameters for EUR population)
2. Allele harmonisation, strand correction and palindromic variant handling
3. Multiple MR estimators (IVW, MR-Egger)
4. Heterogeneity test, horizontal pleiotropy diagnostic
5. External replication with consistent effect-direction filtering
6. TCGA expression exploratory follow-up

All software versions, dataset download links and LD reference panel information are documented within this repository (`sessionInfo.txt`, `data/README_data.txt`).