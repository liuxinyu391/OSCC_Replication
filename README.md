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