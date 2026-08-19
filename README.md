# Multi-tissue Mendelian Randomization Analysis of Oral Cancer Susceptibility

This repository contains the complete analysis pipeline for the manuscript.

## Requirements
- R 4.x (see sessionInfo.txt for full package versions)
- PLINK 1.9 (Windows 64‑bit)

## How to reproduce
1. Clone or download this repository.
2. Edit `WORK_DIR` in `analysis_pipeline.R` to your local path.
3. Download required raw data as described in `data/README_data.txt`.
4. Run `source("analysis_pipeline.R")` in RStudio.

After running the script, local folders `results/` and `results/figures/` will be created automatically to store output tables and figures. These output files are not uploaded to the repository.

## Data availability
All raw data are publicly available (GTEx, FinnGen, GWAS Catalog, TCGA). Due to file size limits, raw data are not included here. See `data/README_data.txt` for download links.

## Results
Output tables and figures will be generated locally under `results/` and `results/figures/` after successful execution of the analysis script. Result files are not hosted in this repository.

## Contact
For questions regarding the code, please contact the corresponding author through the editorial office.

## MR Reporting Standard (STROBE-MR Checklist Compliance)

This analysis follows the STROBE-MR reporting guidelines. All analytical steps required by the checklist are implemented in `analysis_pipeline.R`, including SNP selection, allele harmonisation, multiple MR estimators, heterogeneity and pleiotropy diagnostics, external replication, and TCGA exploratory analysis.

All software versions, dataset download links and LD reference panel information are documented within this repository (`sessionInfo.txt`, `data/README_data.txt`).