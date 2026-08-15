================ Raw data download and placement guide ================

All datasets are publicly available. Please create the subdirectories below
and place the downloaded files accordingly.

1. GTEx v10 eQTLs (four tissues)
   Place these files in data/GTEx_selected/:
   - Whole_Blood.v10.eQTLs.signif_pairs.parquet
   - Muscle_Skeletal.v10.eQTLs.signif_pairs.parquet
   - Esophagus_Mucosa.v10.eQTLs.signif_pairs.parquet
   - Minor_Salivary_Gland.v10.eQTLs.signif_pairs.parquet
   - GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz
   Source: GTEx Portal (https://gtexportal.org/)

2. FinnGen oral cavity cancer GWAS (R12)
   Download from:
   https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_ORALCAVITY_EXALLC.gz
   Place in data/GWAS/ (keep original filename).

3. External validation GWAS
   - GCST90041888 (Pan-UK Biobank)
     File: 34737426-GCST90041888-EFO_0005570.h.tsv.gz
     Download: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90041001-GCST90042000/GCST90041888/harmonised/34737426-GCST90041888-EFO_0005570.h.tsv.gz
   - GCST012237 (Lesseur et al. 2016)
     File: GCST012237_buildGRCh37.tsv
     Download: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST012001-GCST013000/GCST012237/GCST012237_buildGRCh37.tsv
   Place these files in data/external_gwas/ .

4. LD reference panel (1000 Genomes Phase 3, EUR)
   Download from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
   Extract so that data/ld_ref/EUR/EUR.bed, EUR.bim, EUR.fam exist.

5. TCGA-HNSC data
   The script will automatically download expression data via TCGAbiolinks
   into data/GDCdata/. No manual download is needed.

Note: All data are summary-level and publicly available.