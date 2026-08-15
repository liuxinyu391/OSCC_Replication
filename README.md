# Multi-tissue Mendelian Randomization Analysis of Oral Cancer Susceptibility

This repository contains the complete analysis pipeline for the manuscript submitted to *Research Connections*.

## Requirements
- R 4.x (see `sessionInfo.txt` for full package versions)
- PLINK 1.9 (Windows 64-bit)

## How to reproduce
1. Clone or download this repository.
2. Edit `WORK_DIR` in `analysis_pipeline.R` to your local path.
3. Download required raw data as described in `data/README_data.txt`.
4. Run `source("analysis_pipeline.R")` in RStudio.

## Data availability
All raw data are publicly available (GTEx, FinnGen, GWAS Catalog, TCGA). Due to file size limits, raw data are not included here. See `data/README_data.txt` for download links.

## Results
All output files are in `results/` and `results/figures/`.

## Contact
Corresponding author: lizhg@lzu.edu.cn