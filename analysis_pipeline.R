# ============================================================================
# Data Sources and Input Files
# ============================================================================
# The following publicly available datasets are required. Please download and
# place them in the directories indicated before running this script.
#
# 1. GTEx v10 cis-eQTL data (four tissues)
#    Files:
#      Whole_Blood.v10.eQTLs.signif_pairs.parquet
#      Muscle_Skeletal.v10.eQTLs.signif_pairs.parquet
#      Esophagus_Mucosa.v10.eQTLs.signif_pairs.parquet
#      Minor_Salivary_Gland.v10.eQTLs.signif_pairs.parquet
#    Source: GTEx Portal (https://gtexportal.org/)
#    dbGaP accession: phs000424.v10.p1
#    Put in: data/GTEx_selected/
#
# 2. GTEx variant lookup table
#    File:
#      GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz
#    Source: GTEx Portal (https://gtexportal.org/)
#    Put in: data/GTEx_selected/
#
# 3. FinnGen R12 oral cavity cancer GWAS (discovery cohort)
#    File: finngen_R12_C3_ORALCAVITY_EXALLC.gz
#    Phenocode: C3_ORALCAVITY_EXALLC
#    Source: FinnGen (https://www.finngen.fi/)
#    Direct download:
#      https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_ORALCAVITY_EXALLC.gz
#    Put in: data/GWAS/
#
# 4. External validation GWAS
#    (a) GCST90041888 (Pan-UK Biobank; Jiang et al. 2021)
#        GWAS Catalog accession: GCST90041888
#        File: 34737426-GCST90041888-EFO_0005570.h.tsv.gz
#        Direct download:
#          http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90041001-GCST90042000/GCST90041888/harmonised/34737426-GCST90041888-EFO_0005570.h.tsv.gz
#        Put in: data/external_gwas/
#
#    (b) GCST012237 (Lesseur et al. 2016)
#        GWAS Catalog accession: GCST012237
#        File: GCST012237_buildGRCh37.tsv
#        Direct download:
#          https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST012001-GCST013000/GCST012237/GCST012237_buildGRCh37.tsv
#        Put in: data/external_gwas/
#
# 5. LD reference panel (1000 Genomes Phase 3, EUR)
#    File: 1kg.v3.tgz (contains EUR.bed, EUR.bim, EUR.fam)
#    Source: MRC IEU OpenGWAS
#    Direct download: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
#    Put in: data/ld_ref/EUR/
#
# 6. PLINK 1.9 (Windows 64-bit stable version)
#    Source: PLINK (https://www.cog-genomics.org/plink/)
#    The plinkbinr R package will automatically locate it, or you may install
#    it manually and ensure it is on the system PATH.
#
# 7. TCGA-HNSC data
#    Downloaded automatically via TCGAbiolinks, or manually placed in:
#      data/GDCdata/TCGA-HNSC/Transcriptome_Profiling/Gene_Expression_Quantification/
#    Source: GDC Data Portal (https://portal.gdc.cancer.gov/)
#
# Note: All datasets are publicly available summary-level data and contain no
# individual-level or identifiable information.
# ============================================================================

# ========================== 0. Environment Setup ==========================
rm(list = ls())
gc()

# ---------- Set working directory (modify to your actual path)----------
WORK_DIR <- "E:/OSCC_Replication"   # 若使用其他路径，修改此处
setwd(WORK_DIR)

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("results/TCGA", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("data/GTEx_clean", showWarnings = FALSE, recursive = TRUE)
dir.create("data/GWAS", showWarnings = FALSE, recursive = TRUE)
dir.create("data/external_gwas", showWarnings = FALSE, recursive = TRUE)
# ---------- Install and load packages ----------
# 安装 plinkbinr（GitHub 包，单独处理）
if (!require("plinkbinr", character.only = TRUE)) {
  if (!require("devtools", character.only = TRUE)) {
    install.packages("devtools", dependencies = TRUE)
  }
  devtools::install_github("explodecomputer/plinkbinr")
  library(plinkbinr, character.only = TRUE)
}
cran_pkgs <- c("arrow", "data.table", "dplyr", "tidyr", "vroom", "ggplot2", "ggrepel",
               "cowplot", "progress", "boot", "ppcor", "here")
bioc_pkgs <- c("TCGAbiolinks", "DESeq2", "SummarizedExperiment", "edgeR",
               "GSVA", "pheatmap", "survival", "survminer", "clusterProfiler",
               "org.Hs.eg.db", "DOSE", "enrichplot", "BiocManager")
mr_pkgs <- c("TwoSampleMR", "coloc", "plinkbinr", "ieugwasr")

all_pkgs <- c(cran_pkgs, bioc_pkgs, mr_pkgs)
for (pkg in all_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c(bioc_pkgs, "TwoSampleMR", "coloc", "ieugwasr")) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}
cat("✅ 所有包加载完成。All packages loaded.\n")


# ========================== 1. Data Preprocessing  ==========================
cat("\n===== 1. Data Preprocessing =====\n")

# 确保必要目录存在
dir.create("data/GTEx_clean", showWarnings = FALSE, recursive = TRUE)
dir.create("data/GWAS",       showWarnings = FALSE, recursive = TRUE)
dir.create("data/external_gwas", showWarnings = FALSE, recursive = TRUE)


# ---------- 1.1 生成 rsid 查找表的 RDS（只需运行一次，之后秒读）----------
rsid_rds <- "data/GTEx_clean/rsid_lookup.rds"
if (!file.exists(rsid_rds)) {
  cat("正在读取 rsid 查找表并保存为 RDS，之后会快很多...\n")
  rsid <- fread("data/GTEx_selected/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz")
  rsid <- rsid[, c(1:5, 7)]   # 只保留需要的列
  saveRDS(rsid, rsid_rds)
  cat("✅ rsid RDS 已保存。\n")
}

# ---------- 1.2 生成单个组织的 eQTL 合并 CSV ----------


tissues_raw <- list(
  Whole_Blood          = "Whole_Blood",
  Muscle_Skeletal      = "Muscle_Skeletal",
  Esophagus_Mucosa     = "Esophagus_Mucosa",
  Minor_Salivary_Gland = "Minor_Salivary_Gland"
)
# ------ 手动选择要处理的组织Manually sellect the tissues to be processed------
# ★ Permitted values: "Whole_Blood", "Muscle_Skeletal", "Esophagus_Mucosa", "Minor_Salivary_Gland"
SELECTED_TISSUE <- "Esophagus_Mucosa"  

tissue <- SELECTED_TISSUE
out_csv <- file.path("data/GTEx_clean", paste0(tissue, "_all_eqtl.csv"))

if (file.exists(out_csv)) {
  cat(tissue, "的 CSV 已存在，跳过。\n")
} else {
  cat("正在生成", tissue, "的合并 CSV...\n")
  txt_file <- file.path("data/GTEx_selected", 
                        paste0(tissues_raw[[tissue]], ".v10.eQTLs.signif_pairs.txt"))
  
  # 如果 txt 不存在但 parquet 存在，则自动转换（可能较慢，建议提前手动转换）
  if (!file.exists(txt_file)) {
    parquet_file <- file.path("data/GTEx_selected", 
                              paste0(tissues_raw[[tissue]], ".v10.eQTLs.signif_pairs.parquet"))
    if (file.exists(parquet_file)) {
      cat("  未找到 txt，自动从 parquet 转换...（如已手动转好 txt，本步可跳过）\n")
      eqtl_data <- read_parquet(parquet_file)
      fwrite(eqtl_data, txt_file, sep = "\t")
    } else {
      stop("缺少组织 ", tissue, " 的原始数据，请下载并放入 data/GTEx_selected/")
    }
  }
  
  # 读取原始 eQTL 数据
  raw_eqtl <- fread(txt_file)
  
  # 从 RDS 快速读取 rsid 表
  rsid <- readRDS(rsid_rds)
  
  # 合并并清洗
  merged <- merge(rsid, raw_eqtl, by = "variant_id") %>% na.omit() %>% as.data.frame()
  merged$gene_id <- sapply(strsplit(merged$gene_id, "\\."), `[`, 1)
  fwrite(merged, out_csv)
  cat("  ✅", out_csv, "已生成。\n")
}

# ---------- 1.3 将 CSV 转为 RDS（供后续分析使用）----------
prep_eqtl <- function(input_csv, output_rds) {
  if (!file.exists(input_csv)) {
    warning("文件不存在: ", input_csv)
    return(NULL)
  }
  dt <- fread(input_csv)
  if ("rs_id_dbSNP155_GRCh38p13" %in% names(dt))
    setnames(dt, "rs_id_dbSNP155_GRCh38p13", "rsID")
  if ("gene_id" %in% names(dt) && !"gene" %in% names(dt))
    dt[, gene := sapply(strsplit(gene_id, "\\."), `[`, 1)]
  if ("chr" %in% names(dt)) {
    dt[, chr_num := as.integer(gsub("chr", "", chr))]
    dt <- dt[!is.na(chr_num) & chr_num %in% 1:22]
  }
  saveRDS(dt, output_rds)
  cat("✅ 已保存: ", output_rds, "\n")
  return(dt)
}

# eQTL 文件列表
eqtl_files <- list(
  blood   = list(csv = "data/GTEx_clean/Whole_Blood_all_eqtl.csv", 
                 rds = "data/GTEx_clean/Whole_Blood_all_eqtl.rds"),
  eso     = list(csv = "data/GTEx_clean/Esophagus_Mucosa_all_eqtl.csv", 
                 rds = "data/GTEx_clean/Esophagus_Mucosa_all_eqtl.rds"),
  muscle  = list(csv = "data/GTEx_clean/Muscle_Skeletal_all_eqtl.csv", 
                 rds = "data/GTEx_clean/Muscle_Skeletal_all_eqtl.rds"),
  sal     = list(csv = "data/GTEx_clean/Minor_Salivary_Gland_all_eqtl.csv", 
                 rds = "data/GTEx_clean/Minor_Salivary_Gland_all_eqtl.rds")
)
for (f in eqtl_files) {
  if (file.exists(f$csv)) prep_eqtl(f$csv, f$rds)
}

# ---------- 1.4 其他数据集处理 ----------
# FinnGen GWAS
if (file.exists("data/GWAS/finngen_R12_C3_ORALCAVITY_EXALLC.gz")) {
  outcome <- vroom("data/GWAS/finngen_R12_C3_ORALCAVITY_EXALLC.gz", show_col_types = FALSE) %>% as.data.table()
  if ("#chrom" %in% names(outcome)) setnames(outcome, "#chrom", "chrom")
  outcome[, chrom_num := as.integer(gsub("chr", "", chrom))]
  saveRDS(outcome, "data/GWAS/finngen_R12_clean.rds")
  cat("✅ FinnGen GWAS 已保存。\n")
} else {
  stop("❌ 缺少 FinnGen GWAS 文件，请下载后放置于 data/GWAS/")
}

# 外部 GWAS
prep_ext_gwas <- function(input_path, output_rds) {
  if (!file.exists(input_path)) { warning("外部GWAS不存在: ", input_path); return(NULL) }
  dt <- fread(input_path,
              select = c("variant_id","beta","standard_error","p_value","effect_allele","other_allele","effect_allele_frequency"),
              quote = "",
              fill = TRUE)
  numeric_cols <- c("beta","standard_error","p_value","effect_allele_frequency")
  dt[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]
  
  setnames(dt, old = names(dt), new = c("SNP","beta","se","pval","effect_allele","other_allele","eaf"))
  saveRDS(dt, output_rds)
  cat("✅ 已保存: ", output_rds, "\n")
  return(dt)
}
prep_ext_gwas("data/external_gwas/34737426-GCST90041888-EFO_0005570.h.tsv.gz", "data/external_gwas/GCST90041888.rds")
prep_ext_gwas("data/external_gwas/GCST012237_buildGRCh37.tsv", "data/external_gwas/GCST012237.rds")

# 加载所有预处理数据（后续分析直接使用）
outcome <- readRDS("data/GWAS/finngen_R12_clean.rds")
eqtl_blood <- readRDS("data/GTEx_clean/Whole_Blood_all_eqtl.rds")
eqtl_eso   <- readRDS("data/GTEx_clean/Esophagus_Mucosa_all_eqtl.rds")
eqtl_muscle<- readRDS("data/GTEx_clean/Muscle_Skeletal_all_eqtl.rds")
eqtl_sal   <- readRDS("data/GTEx_clean/Minor_Salivary_Gland_all_eqtl.rds")
ext_gwas1  <- readRDS("data/external_gwas/GCST90041888.rds")
ext_gwas2  <- readRDS("data/external_gwas/GCST012237.rds")

plink_exe <- tryCatch(get_plink_exe(), error = function(e) NULL)
if (is.null(plink_exe) || Sys.which(plink_exe) == "") stop("❌ PLINK 未找到。")
ld_ref_dir <- "data/ld_ref/EUR/EUR"
if (!file.exists(paste0(ld_ref_dir, ".bed"))) stop("❌ LD 参考文件不存在: ", ld_ref_dir, ".bed")

cat("✅ 数据预处理完成。\n")
# ========================== 2. Core MR Analysis (Batch-wise with Checkpoint/Resume） ==========================
cat("\n===== 2. 核心 MR 分析 =====\n")

run_mr_tissue <- function(org_name, eqtl_data, pval_thresh, ld_kb = 1000, ld_r2 = 0.01,
                          batch_size = 50) {
  cat("\n组织:", org_name, "\n")
  eqtl_sub <- eqtl_data[pval_nominal < pval_thresh, ]
  if (nrow(eqtl_sub) == 0) { cat("❌ 无显著 eQTL。\n"); return(NULL) }
  gene_list <- unique(eqtl_sub$gene)
  cat("Total genes:", length(gene_list), "\n")
  
  output_file <- file.path("results", paste0("MR_results_", org_name, ".csv"))
  progress_file <- file.path("results", paste0("MR_progress_", org_name, ".txt"))
  
  # 清理无效进度
  if (file.exists(progress_file) && !file.exists(output_file)) {
    cat("⚠️ 检测到进度文件但无结果文件，将重新开始分析。\n")
    file.remove(progress_file)
  }
  
  done_genes <- c()
  if (file.exists(progress_file)) {
    done_genes <- readLines(progress_file)
    cat("Skipped", length(done_genes), "completed genes\n")
  }
  
  remaining_genes <- setdiff(gene_list, done_genes)
  if (length(remaining_genes) == 0) {
    cat("所有基因已完成。\n")
    if (file.exists(output_file)) return(fread(output_file)) else return(NULL)
  }
  
  n_batches <- ceiling(length(remaining_genes) / batch_size)
  cat("Will be divided into", n_batches, "batches, up to", batch_size, "genes per batch\n")
  
  # ---------- 内嵌单基因 MR 函数（含加权中位数和 MR-Egger） ----------
  mr_one_gene <- function(exp_sub) {
    g <- unique(exp_sub$gene)[1]
    if (nrow(exp_sub) == 0) return(NULL)
    exp_fmt <- tryCatch(format_data(as.data.frame(exp_sub), type="exposure", phenotype_col="gene",
                                    snp_col="rsID", beta_col="slope", se_col="slope_se", eaf_col="af",
                                    effect_allele_col="alt", other_allele_col="ref", pval_col="pval_nominal",
                                    chr_col="chr", pos_col="pos"), error=function(e) NULL)
    if (is.null(exp_fmt) || nrow(exp_fmt) < 1) return(NULL)
    exp_fmt$id.exposure <- g
    clump_input <- data.frame(rsid=exp_fmt$SNP, pval=exp_fmt$pval.exposure)
    clumped <- tryCatch(ieugwasr::ld_clump(clump_input, clump_kb=ld_kb, clump_r2=ld_r2,
                                           clump_p=1, bfile=ld_ref_dir, plink_bin=plink_exe, pop="EUR"),
                        error=function(e) NULL)
    if (is.null(clumped) || nrow(clumped)==0) return(NULL)
    exp_clump <- exp_fmt[exp_fmt$SNP %in% clumped$rsid, ]; if (nrow(exp_clump)<1) return(NULL)
    out_raw <- outcome[rsids %in% exp_clump$SNP, ]; if (nrow(out_raw)==0) return(NULL)
    out_fmt <- format_data(as.data.frame(out_raw), type="outcome", snp_col="rsids", beta_col="beta",
                           se_col="sebeta", effect_allele_col="alt", other_allele_col="ref",
                           eaf_col="af_alt", pval_col="pval")
    out_fmt$id.outcome <- "OSCC"
    mr_data <- harmonise_data(exp_clump, out_fmt); if (nrow(mr_data)==0) return(NULL)
    mr_data$F <- (mr_data$beta.exposure / mr_data$se.exposure)^2; minF <- min(mr_data$F, na.rm=TRUE)
    
    # 请求所有 MR 方法（IVW、Wald ratio、加权中位数、MR-Egger）
    mr_res <- mr(mr_data, method_list=c("mr_ivw","mr_wald_ratio",
                                        "mr_weighted_median","mr_egger_regression"))
    
    # 提取 IVW（或 Wald ratio 作为单 SNP 的替代）
    ivw <- mr_res[mr_res$method=="Inverse variance weighted",]
    if (nrow(ivw)==0) { ivw <- mr_res[mr_res$method=="Wald ratio",]; if (nrow(ivw)==0) return(NULL) }
    
    # 提取加权中位数（若存在）
    wm <- mr_res[mr_res$method=="Weighted median",]
    wm_b <- if (nrow(wm)>0) wm$b[1] else NA
    wm_se <- if (nrow(wm)>0) wm$se[1] else NA
    wm_p <- if (nrow(wm)>0) wm$pval[1] else NA
    
    # 提取 MR-Egger（若存在）
    egger <- mr_res[mr_res$method=="MR Egger",]
    egger_b <- if (nrow(egger)>0) egger$b[1] else NA
    egger_se <- if (nrow(egger)>0) egger$se[1] else NA
    egger_p <- if (nrow(egger)>0) egger$pval[1] else NA
    
    # 异质性
    het <- mr_heterogeneity(mr_data)
    het_p <- ifelse(nrow(het)>0, het$Q_pval[1], NA)
    
    # 多效性（MR-Egger 截距）
    pleio_p <- NA
    if (nrow(mr_data)>=3) {
      pleio <- tryCatch(mr_pleiotropy_test(mr_data), error=function(e) NULL)
      pleio_p <- if (!is.null(pleio)) pleio$pval else NA
    }
    
    # 返回结果（保留原有列，并新增敏感性分析列）
    data.table(
      gene = g,
      tissue = org_name,
      nsnp = nrow(mr_data),
      minF = minF,
      ivw_b = ivw$b,
      ivw_se = ivw$se,
      ivw_pval = ivw$pval,
      ivw_OR = exp(ivw$b),
      ivw_OR_lci = exp(ivw$b - 1.96*ivw$se),
      ivw_OR_uci = exp(ivw$b + 1.96*ivw$se),
      weighted_median_b = wm_b,
      weighted_median_se = wm_se,
      weighted_median_pval = wm_p,
      egger_b = egger_b,
      egger_se = egger_se,
      egger_pval = egger_p,
      het_p = het_p,
      pleio_p = pleio_p
    )
  }
  # ------------------------------------------------------------
  
  for (batch_idx in 1:n_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, length(remaining_genes))
    batch_genes <- remaining_genes[start_idx:end_idx]
    
    cat(sprintf("\n===== 批次 %d/%d (基因 %d-%d) =====\n", 
                batch_idx, n_batches, start_idx, end_idx))
    
    batch_sub <- eqtl_sub[gene %in% batch_genes, ]
    genes_split <- split(batch_sub, by = "gene")
    
    res_list <- list()
    success_genes <- c()
    for (i in seq_along(genes_split)) {
      g <- names(genes_split)[i]
      res <- tryCatch(mr_one_gene(genes_split[[i]]), error = function(e) NULL)
      if (!is.null(res)) {
        res_list[[length(res_list) + 1]] <- res
        success_genes <- c(success_genes, g)
      }
    }
    
    if (length(res_list) > 0) {
      batch_df <- rbindlist(res_list, fill = TRUE)
      fwrite(batch_df, output_file, append = file.exists(output_file))
    }
    
    done_genes <- unique(c(done_genes, success_genes))
    writeLines(done_genes, progress_file)
    
    if (length(success_genes) < length(batch_genes)) {
      cat(sprintf("  本批 %d 个基因中 %d 个成功，%d 个将在下次重试\n",
                  length(batch_genes), length(success_genes), length(batch_genes) - length(success_genes)))
    }
    
    rm(res_list, batch_df, batch_sub, genes_split)
    gc()
    Sys.sleep(1)
  }
  
  if (file.exists(output_file)) {
    cat("✅", org_name, "完成，结果已保存。\n")
    return(fread(output_file))
  } else {
    cat("❌", org_name, "未产生任何有效 MR 结果。\n")
    return(NULL)
  }
}

# ---------- Call (Tier 1 strict parameters, Tier 2 lenient parameters）----------
mr_blood  <- run_mr_tissue("Whole_Blood",         eqtl_blood,  5e-8, ld_kb=1000, ld_r2=0.01, batch_size=50)
mr_eso    <- run_mr_tissue("Esophagus_Mucosa",     eqtl_eso,    5e-8, ld_kb=1000, ld_r2=0.01, batch_size=50)
mr_muscle <- run_mr_tissue("Muscle_Skeletal",      eqtl_muscle, 5e-8, ld_kb=1000, ld_r2=0.01, batch_size=50)
mr_sal    <- run_mr_tissue("Minor_Salivary_Gland", eqtl_sal,    1e-6, ld_kb=500,  ld_r2=0.1,  batch_size=100)

# Safely read empty results
read_or_empty <- function(file) if (file.exists(file)) fread(file) else data.table()
mr_blood  <- if (is.null(mr_blood))  read_or_empty("results/MR_results_Whole_Blood.csv") else mr_blood
mr_eso    <- if (is.null(mr_eso))    read_or_empty("results/MR_results_Esophagus_Mucosa.csv") else mr_eso
mr_muscle <- if (is.null(mr_muscle)) read_or_empty("results/MR_results_Muscle_Skeletal.csv") else mr_muscle
mr_sal    <- if (is.null(mr_sal))    read_or_empty("results/MR_results_Minor_Salivary_Gland.csv") else mr_sal

library(data.table)

# Read MR results for four tissues
mr_blood  <- fread("results/MR_results_Whole_Blood.csv")
mr_eso    <- fread("results/MR_results_Esophagus_Mucosa.csv")
mr_muscle <- fread("results/MR_results_Muscle_Skeletal.csv")
mr_sal    <- fread("results/MR_results_Minor_Salivary_Gland.csv")

# ========================== 3. Evidence Integration ==========================
cat("\n===== 3. 证据整合 =====\n")
integrate_evidence <- function(mr_blood, mr_eso, mr_muscle, mr_sal, pval_type="fdr", suffix="") {
  main_list <- list(Whole_Blood=mr_blood, Esophagus_Mucosa=mr_eso, Muscle_Skeletal=mr_muscle)
  for (tissue in names(main_list)) {
    df <- main_list[[tissue]]
    if (pval_type=="fdr") { df[, ivw_fdr:=p.adjust(ivw_pval, method="fdr")]; df[, strong:=(ivw_fdr<0.05 & minF>10 & (is.na(het_p)|het_p>0.05) & (is.na(pleio_p)|pleio_p>0.05))] }
    else { df[, strong:=(ivw_pval<0.05 & minF>10 & (is.na(het_p)|het_p>0.05) & (is.na(pleio_p)|pleio_p>0.05))] }
    main_list[[tissue]] <- df
  }
  sal_support <- mr_sal[, .(gene, sal_b=ivw_b, sal_support=ivw_pval<0.05)]
  all_genes <- unique(c(main_list$Whole_Blood$gene, main_list$Esophagus_Mucosa$gene, main_list$Muscle_Skeletal$gene, sal_support$gene))
  integrated <- data.table(gene=all_genes)
  eso_df <- main_list$Esophagus_Mucosa[, .(gene, eso_strong=strong, eso_b=ivw_b)]
  integrated <- merge(integrated, eso_df, by="gene", all.x=TRUE)
  integrated <- merge(integrated, sal_support, by="gene", all.x=TRUE)
  integrated[, direction_consistent:=sign(eso_b)==sign(sal_b) & !is.na(eso_b) & !is.na(sal_b)]
  integrated[, sal_final:=sal_support & direction_consistent]
  integrated[, grade:=fcase(eso_strong==TRUE & sal_final==TRUE, "A", eso_strong==TRUE & (is.na(sal_final)|sal_final==FALSE), "B", (is.na(eso_strong)|eso_strong==FALSE) & sal_final==TRUE, "C", default="Exclude")]
  blood_strong <- main_list$Whole_Blood[strong==TRUE, gene]; muscle_strong <- main_list$Muscle_Skeletal[strong==TRUE, gene]
  integrated[grade=="Exclude" & gene %in% union(blood_strong, muscle_strong), grade:="Exclude (blood/muscle only)"]
  fwrite(integrated, file.path("results", paste0("integrated_genes_", pval_type, suffix, ".csv")))
  return(integrated)
}

integrated_main <- integrate_evidence(mr_blood, mr_eso, mr_muscle, mr_sal, pval_type="fdr", suffix="_main")
integrated_exploratory <- integrate_evidence(mr_blood, mr_eso, mr_muscle, mr_sal, pval_type="nominal", suffix="_exploratory")
integrated <- fread("results/integrated_genes_nominal_exploratory.csv")
AB_genes <- integrated[grade %in% c("A","B"), gene]; C_genes <- integrated[grade=="C", gene]

library(org.Hs.eg.db)
library(clusterProfiler)
library(data.table)

all_genes_ens <- unique(integrated$gene)
genes_clean <- sub("\\..*", "", all_genes_ens)
anno <- bitr(genes_clean, fromType = "ENSEMBL",
             toType = c("SYMBOL", "ENTREZID", "GENETYPE"),
             OrgDb = org.Hs.eg.db)
setDT(anno)
setnames(anno, "ENSEMBL", "gene")
fwrite(anno, "results/gene_annotation.csv")
cat("✅ 基因注释文件已保存至 results/gene_annotation.csv\n")

# ========================== 4. Co-localization analysis============================================== 
cat("\n===== 4. Co-localization analysis =====\n")

# ---------- 4.1 Single-gene colocalization function ----------
coloc_one_gene <- function(g, eqtl_data, sample_size) {
  gene_eqtl <- eqtl_data[gene == g, ]
  if (nrow(gene_eqtl) == 0) return(NULL)
  chr_counts <- table(gene_eqtl$chr_num)
  if (length(chr_counts) == 0) return(NULL)
  main_chr <- as.integer(names(chr_counts)[which.max(chr_counts)])
  gene_eqtl <- gene_eqtl[chr_num == main_chr]
  if (nrow(gene_eqtl) == 0) return(NULL)
  min_pos <- min(gene_eqtl$pos, na.rm = TRUE)
  max_pos <- max(gene_eqtl$pos, na.rm = TRUE)
  region_start <- max(0, min_pos - 500000)
  region_end <- max_pos + 500000
  
  eqtl_region <- eqtl_data[chr_num == main_chr & pos >= region_start & pos <= region_end,
                           .(rsID, pos, pval_nominal, af)]
  gwas_region <- outcome[chrom_num == main_chr & pos >= region_start & pos <= region_end,
                         .(rsids, pos, pval, af_alt)]
  if (nrow(eqtl_region) == 0 || nrow(gwas_region) == 0) return(NULL)
  setkey(eqtl_region, pos)
  setkey(gwas_region, pos)
  merged <- merge(eqtl_region, gwas_region, by = "pos", allow.cartesian = FALSE)
  merged <- merged[rsID == rsids & !is.na(af) & !is.na(af_alt)]
  merged <- merged[!duplicated(rsID)]
  if (nrow(merged) < 10) return(NULL)
  
  coloc_res <- tryCatch(
    coloc.abf(
      dataset1 = list(
        snp = merged$rsID, pvalues = merged$pval_nominal,
        N = sample_size, MAF = merged$af, type = "quant"
      ),
      dataset2 = list(
        snp = merged$rsID, pvalues = merged$pval,
        N = 1614 + 378749, MAF = merged$af_alt,
        type = "cc", s = 1614 / (1614 + 378749)
      )
    ),
    error = function(e) NULL
  )
  if (is.null(coloc_res)) return(NULL)
  data.table(
    gene = g,
    nsnp_coloc = as.numeric(coloc_res$summary["nsnps"]),
    PP.H4 = as.numeric(coloc_res$summary["PP.H4.abf"])
  )
}

# ---------- 4.2 Batch colocalization funcation ----------
run_coloc_batch <- function(gene_list, eqtl_data, sample_size, label, batch_size = 30) {
  if (length(gene_list) == 0) return(NULL)
  output_file <- file.path("results", paste0("coloc_", label, ".csv"))
  progress_file <- file.path("results", paste0("coloc_progress_", label, ".txt"))
  
  if (file.exists(progress_file) && !file.exists(output_file)) {
    cat("⚠️ 清理无效进度\n")
    file.remove(progress_file)
  }
  done_genes <- if (file.exists(progress_file)) readLines(progress_file) else character(0)
  remaining_genes <- setdiff(gene_list, done_genes)
  if (length(remaining_genes) == 0) {
    if (file.exists(output_file)) return(fread(output_file)) else return(NULL)
  }
  
  n_batches <- ceiling(length(remaining_genes) / batch_size)
  for (batch_idx in 1:n_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, length(remaining_genes))
    batch_genes <- remaining_genes[start_idx:end_idx]
    
    res_list <- list()
    success_genes <- c()
    for (g in batch_genes) {
      res <- tryCatch(coloc_one_gene(g, eqtl_data, sample_size), error = function(e) NULL)
      if (!is.null(res)) {
        res_list[[length(res_list) + 1]] <- res
        success_genes <- c(success_genes, g)
      }
      gc()
    }
    if (length(res_list) > 0) {
      batch_df <- rbindlist(res_list, fill = TRUE)
      fwrite(batch_df, output_file, append = file.exists(output_file))
      rm(batch_df)
    }
    done_genes <- unique(c(done_genes, success_genes))
    writeLines(done_genes, progress_file)
    rm(res_list)
    gc()
    Sys.sleep(1)
  }
  if (file.exists(output_file)) return(fread(output_file)) else return(NULL)
}

# ---------- 4.3 Prepare each gene set ----------
# 排除集基因
excluded_genes <- integrated[grade == "Exclude (blood/muscle only)", gene]
blood_sig <- mr_blood[ivw_pval < 0.05, unique(gene)]
muscle_sig <- mr_muscle[ivw_pval < 0.05, unique(gene)]
excl_blood <- intersect(excluded_genes, blood_sig)
excl_muscle <- intersect(excluded_genes, muscle_sig)

# FDR-C 级基因（180个，来自FDR整合文件）
if (file.exists("results/integrated_genes_fdr_main.csv")) {
  integrated_fdr <- fread("results/integrated_genes_fdr_main.csv")
  C_fdr_genes <- integrated_fdr[grade == "C", gene]
} else {
  stop("缺少 FDR 整合文件 integrated_genes_fdr_main.csv，请检查。")
}

# 名义-C 级基因（63个，来自名义显著整合表，作为对照）
C_nominal <- integrated[grade == "C", gene]

# ---------- 4.4 Run gene-level colocalization for each gene ----------
cat("开始共定位分析...\n")
coloc_AB            <- run_coloc_batch(AB_genes, eqtl_eso, 500, "AB", batch_size = 20)
coloc_blood_excl    <- run_coloc_batch(excl_blood, eqtl_blood, 670, "blood_excluded", batch_size = 30)
coloc_muscle_excl   <- run_coloc_batch(excl_muscle, eqtl_muscle, 500, "muscle_excluded", batch_size = 30)
coloc_C_fdr         <- run_coloc_batch(C_fdr_genes, eqtl_sal, 100, "C_fdr", batch_size = 30)
coloc_C_nominal     <- run_coloc_batch(C_nominal, eqtl_sal, 100, "C_nominal", batch_size = 30)

# ---------- 4.5 Summarize statistical table----------
summarize_set <- function(coloc_df, total_genes, set_name) {
  coloc_success <- if (!is.null(coloc_df)) nrow(coloc_df) else 0
  coloc_sup    <- if (!is.null(coloc_df)) sum(coloc_df$PP.H4 >= 0.5, na.rm = TRUE) else 0
  data.table(
    GeneSet      = set_name,
    TotalGenes   = total_genes,
    ColocSuccess = coloc_success,
    ColocSupport = coloc_sup
  )
}

table3 <- rbind(
  summarize_set(coloc_AB,          length(AB_genes),        "A级+B级（食管黏膜）"),
  summarize_set(coloc_blood_excl,  length(excl_blood),      "血液排除集"),
  summarize_set(coloc_muscle_excl, length(excl_muscle),     "肌肉排除集"),
  summarize_set(coloc_C_fdr,       length(C_fdr_genes),     "FDR-C级（唾液腺，n=180）"),
  summarize_set(coloc_C_nominal,   length(C_nominal),       "名义-C级（唾液腺，n=63）")
)

fwrite(table3, "results/Table3_Summary.csv")
print(table3)
cat("✅ 共定位分析及汇总表格全部完成\n")

if (!exists("anno") || is.null(anno) || !is.data.table(anno)) {
  if (file.exists("results/gene_annotation.csv")) {
    anno <- fread("results/gene_annotation.csv")
    setDT(anno)
    cat("✅ 已加载注释数据\n")
  } else {
    stop("❌ 请先运行生成 anno 的代码，或确保 results/gene_annotation.csv 存在。")
  }
}
if (!all(c("gene", "SYMBOL") %in% names(anno))) {
  stop("❌ anno 缺少 'gene' 或 'SYMBOL' 列")
}

#  1. A级基因小唾液腺共定位 
cat("\n=== 补充：A级基因小唾液腺共定位 ===\n")
A_genes <- integrated[grade == "A", gene]
coloc_A_sal <- NULL
if (length(A_genes) > 0) {
  coloc_A_sal <- run_coloc_batch(A_genes, eqtl_sal, 100, "A_sal")
  if (!is.null(coloc_A_sal) && nrow(coloc_A_sal) > 0) {
    cols_keep <- setdiff(names(coloc_A_sal), grep("^SYMBOL", names(coloc_A_sal), value = TRUE))
    coloc_A_sal_clean <- coloc_A_sal[, ..cols_keep]
    anno_sub <- anno[, .(gene, SYMBOL)]
    coloc_A_sal_anno <- merge(coloc_A_sal_clean, anno_sub, by = "gene", all.x = TRUE)
    # 去重
    coloc_A_sal_anno <- coloc_A_sal_anno[, .SD[which.max(PP.H4)], by = gene]
    fwrite(coloc_A_sal_anno, "results/coloc_A_sal.csv")
    cat("✅ A级唾液腺共定位结果已保存: results/coloc_A_sal.csv\n")
  } else {
    cat("⚠️ A级基因共定位无结果\n")
  }
} else {
  cat("⚠️ 无A级基因\n")
}

#  2. B级基因小唾液腺共定位 
cat("\n=== 补充：B级基因小唾液腺共定位 ===\n")
B_genes <- integrated[grade == "B", gene]
coloc_B_sal <- NULL
if (length(B_genes) > 0) {
  coloc_B_sal <- run_coloc_batch(B_genes, eqtl_sal, 100, "B_sal")
  if (!is.null(coloc_B_sal) && nrow(coloc_B_sal) > 0) {
    cols_keep <- setdiff(names(coloc_B_sal), grep("^SYMBOL", names(coloc_B_sal), value = TRUE))
    coloc_B_sal_clean <- coloc_B_sal[, ..cols_keep]
    anno_sub <- anno[, .(gene, SYMBOL)]
    coloc_B_sal_anno <- merge(coloc_B_sal_clean, anno_sub, by = "gene", all.x = TRUE)
    # 去重
    coloc_B_sal_anno <- coloc_B_sal_anno[, .SD[which.max(PP.H4)], by = gene]
    fwrite(coloc_B_sal_anno, "results/coloc_B_sal.csv")
    cat("✅ B级唾液腺共定位结果已保存: results/coloc_B_sal.csv\n")
  } else {
    cat("⚠️ B级基因共定位无结果\n")
  }
} else {
  cat("⚠️ 无B级基因\n")
}

#  3. A vs B 共定位支持率 Fisher 检验 
cat("\n=== A vs B 共定位支持率比较 ===\n")
# 重新读取去重后的文件（确保数据正确）
if (file.exists("results/coloc_A_sal.csv")) {
  coloc_A_sal <- fread("results/coloc_A_sal.csv")
  setDT(coloc_A_sal)
}
if (file.exists("results/coloc_B_sal.csv")) {
  coloc_B_sal <- fread("results/coloc_B_sal.csv")
  setDT(coloc_B_sal)
}

# 强制过滤，只保留属于当前分级的基因
coloc_A_valid <- coloc_A_sal[gene %in% A_genes, ]
coloc_B_valid <- coloc_B_sal[gene %in% B_genes, ]

A_total <- nrow(coloc_A_valid)
B_total <- nrow(coloc_B_valid)
A_support <- sum(coloc_A_valid$PP.H4 >= 0.5, na.rm = TRUE)
B_support <- sum(coloc_B_valid$PP.H4 >= 0.5, na.rm = TRUE)

cat("A级基因总数:", length(A_genes), "\n")
cat("A级共定位成功数:", A_total, "\n")
cat("B级基因总数:", length(B_genes), "\n")
cat("B级共定位成功数:", B_total, "\n")

if (A_total > 0 && B_total > 0) {
  mat <- matrix(c(A_support, A_total - A_support,
                  B_support, B_total - B_support), nrow = 2, byrow = TRUE)
  fisher_res <- fisher.test(mat)
  cat(sprintf("A级支持率: %d/%d = %.1f%%\n", A_support, A_total, 100*A_support/A_total))
  cat(sprintf("B级支持率: %d/%d = %.1f%%\n", B_support, B_total, 100*B_support/B_total))
  cat("Fisher 检验 P =", fisher_res$p.value, "\n")
  fisher_out <- data.table(Comparison = "A_vs_B_Support", 
                           A_rate = paste0(A_support, "/", A_total),
                           B_rate = paste0(B_support, "/", B_total),
                           P_value = fisher_res$p.value)
  fwrite(fisher_out, "results/A_vs_B_coloc_fisher.csv")
} else {
  cat("⚠️ A或B级共定位数据不足，无法比较\n")
}

#  4. 更新表3（添加A/B级共定位统计行） 
cat("\n=== 更新表3：添加A/B级共定位统计 ===\n")
if (file.exists("results/Table3_Summary.csv")) {
  table3 <- fread("results/Table3_Summary.csv")
} else {
  table3 <- data.table()
}

new_rows <- rbind(
  data.table(GeneSet = "A级（唾液腺）", 
             TotalGenes = length(A_genes),
             ColocSuccess = A_total,
             ColocSupport = A_support,
             StrongSupport = sum(coloc_A_valid$PP.H4 > 0.8, na.rm = TRUE)),
  data.table(GeneSet = "B级（唾液腺）",
             TotalGenes = length(B_genes),
             ColocSuccess = B_total,
             ColocSupport = B_support,
             StrongSupport = sum(coloc_B_valid$PP.H4 > 0.8, na.rm = TRUE))
)

table3 <- table3[!GeneSet %in% c("A级（唾液腺）", "B级（唾液腺）")]
table3 <- rbind(table3, new_rows, fill = TRUE)
fwrite(table3, "results/Table3_Summary_updated.csv")
cat("✅ 更新后的表3已保存: results/Table3_Summary_updated.csv\n")
print(table3)

#  5. B级基因选择偏倚分析
cat("\n=== B级基因选择偏倚分析（补充表15a）最终修正版 ===\n")

B_genes_all <- integrated[grade == "B", gene]
success_genes <- unique(coloc_B_valid$gene)
fail_genes <- setdiff(B_genes_all, success_genes)
cat("B级成功分析基因数:", length(success_genes), "\n")
cat("B级失败分析基因数:", length(fail_genes), "\n")

if (length(success_genes) > 0 || length(fail_genes) > 0) {
  # 基因位置信息
  if (!exists("eqtl_sal") || is.null(eqtl_sal)) {
    eqtl_sal <- readRDS("data/GTEx_clean/Minor_Salivary_Gland_all_eqtl.rds")
  }
  eqtl_sal[, chr_num := as.integer(gsub("chr", "", chr))]
  eqtl_sal <- eqtl_sal[!is.na(chr_num) & chr_num %in% 1:22]
  eqtl_sal[, gene_clean := sub("\\..*", "", gene)]
  gene_pos <- eqtl_sal[, .(chromosome = first(chr_num),
                           min_pos = min(pos, na.rm = TRUE),
                           max_pos = max(pos, na.rm = TRUE)), by = gene_clean]
  setnames(gene_pos, "gene_clean", "gene")
  gene_pos[, gene_length := max_pos - min_pos]
  gene_pos[, is_mhc := ifelse(chromosome == 6 & min_pos >= 28500000 & max_pos <= 33500000, TRUE, FALSE)]
  
  # 合并分组
  all_genes <- data.table(gene = c(success_genes, fail_genes),
                          group = c(rep("Success", length(success_genes)), 
                                    rep("Fail", length(fail_genes))))
  all_genes[, gene_clean := sub("\\..*", "", gene)]
  gene_info <- merge(all_genes, gene_pos, by.x = "gene_clean", by.y = "gene", all.x = TRUE)
  gene_info_clean <- gene_info[!is.na(chromosome)]
  
  # 合并 GENETYPE（关键：使用正确的 "protein-coding" 字符串）
  anno[, gene_clean := sub("\\..*", "", gene)]
  gene_info_clean <- merge(gene_info_clean, anno[, .(gene_clean, GENETYPE)], 
                           by = "gene_clean", all.x = TRUE)
  
  # 使用正确的蛋白编码字符串 "protein-coding"（带连字符）
  protein_coding_string <- "protein-coding"
  pc_success <- sum(gene_info_clean[group == "Success", GENETYPE == protein_coding_string], na.rm = TRUE)
  pc_fail <- sum(gene_info_clean[group == "Fail", GENETYPE == protein_coding_string], na.rm = TRUE)
  total_success <- length(success_genes)
  total_fail <- length(fail_genes)
  pc_ratio_success <- pc_success / total_success * 100
  pc_ratio_fail <- pc_fail / total_fail * 100
  
  cat("成功组蛋白编码:", pc_success, "/", total_success, "=", round(pc_ratio_success, 1), "%\n")
  cat("失败组蛋白编码:", pc_fail, "/", total_fail, "=", round(pc_ratio_fail, 1), "%\n")
  
  pc_pval <- fisher.test(matrix(c(pc_success, total_success - pc_success,
                                  pc_fail, total_fail - pc_fail), nrow = 2))$p.value
  
  # 基因长度检验
  gene_len_clean <- gene_info_clean[!is.na(gene_length) & gene_length > 0]
  if (nrow(gene_len_clean) > 0) {
    wilcox_p <- wilcox.test(gene_length ~ group, data = gene_len_clean)$p.value
    median_success <- median(gene_len_clean[group == "Success", gene_length] / 1000, na.rm = TRUE)
    median_fail <- median(gene_len_clean[group == "Fail", gene_length] / 1000, na.rm = TRUE)
    iqr_success <- IQR(gene_len_clean[group == "Success", gene_length] / 1000, na.rm = TRUE)
    iqr_fail <- IQR(gene_len_clean[group == "Fail", gene_length] / 1000, na.rm = TRUE)
  } else {
    wilcox_p <- median_success <- median_fail <- iqr_success <- iqr_fail <- NA
  }
  
  # MHC区域比例
  mhc_success <- sum(gene_info_clean[group == "Success", is_mhc], na.rm = TRUE)
  mhc_fail <- sum(gene_info_clean[group == "Fail", is_mhc], na.rm = TRUE)
  mhc_pval <- fisher.test(matrix(c(mhc_success, total_success - mhc_success,
                                   mhc_fail, total_fail - mhc_fail), nrow = 2))$p.value
  
  # 染色体分布检验
  chr_count <- gene_info_clean[, .N, by = .(group, chromosome)]
  chr_wide <- dcast(chr_count, chromosome ~ group, value.var = "N", fill = 0)
  chisq_p <- tryCatch(chisq.test(chr_wide[, .(Success, Fail)], simulate.p.value = TRUE)$p.value,
                      error = function(e) NA)
  
  # 生成补充表15a
  table15a <- data.frame(
    Feature = c("染色体分布 (χ²检验)",
                "基因长度 (kb, 中位数 [IQR])",
                "位于MHC区域比例 (%)",
                "蛋白编码基因比例 (%)"),
    Success_group = c(if (is.na(chisq_p)) "NA" else paste0("P = ", round(chisq_p, 4)),
                      paste0(round(median_success, 1), " [", round(iqr_success, 1), "]"),
                      paste0(round(100*mhc_success/total_success, 1), "% (", mhc_success, "/", total_success, ")"),
                      paste0(round(pc_ratio_success, 1), "% (", pc_success, "/", total_success, ")")),
    Fail_group = c("",
                   paste0(round(median_fail, 1), " [", round(iqr_fail, 1), "]"),
                   paste0(round(100*mhc_fail/total_fail, 1), "% (", mhc_fail, "/", total_fail, ")"),
                   paste0(round(pc_ratio_fail, 1), "% (", pc_fail, "/", total_fail, ")")),
    P_value = c(round(chisq_p, 4), round(wilcox_p, 4), round(mhc_pval, 4), round(pc_pval, 4))
  )
  dir.create("results/selection_bias_B", showWarnings = FALSE)
  fwrite(table15a, "results/selection_bias_B/Supplementary_Table15a_B_selection_bias.csv")
  cat("✅ 补充表15a已保存至 results/selection_bias_B/Supplementary_Table15a_B_selection_bias.csv\n")
  print(table15a)
} else {
  cat("⚠️ 无B级基因或共定位数据，跳过偏倚分析\n")
}

cat("\n========== 所有补充分析完成 ==========\n")
# ========================== 5. External validation ==========================
# ---------- 5.1 Single-gene validation function ----------
validate_gene <- function(g, eqtl_data, gwas_data) {
  exp_sub <- eqtl_data[gene == g, ]
  if (nrow(exp_sub) == 0) return(NULL)
  exp_fmt <- format_data(
    as.data.frame(exp_sub),
    type = "exposure", phenotype_col = "gene",
    snp_col = "rsID", beta_col = "slope", se_col = "slope_se", eaf_col = "af",
    effect_allele_col = "alt", other_allele_col = "ref", pval_col = "pval_nominal",
    chr_col = "chr", pos_col = "pos"
  )
  exp_fmt$id.exposure <- g
  if (nrow(exp_fmt) < 1) return(NULL)
  
  # LD clumping（与 Tier 1 严格参数一致：1000 kb, r² < 0.01）
  clump_input <- data.frame(rsid = exp_fmt$SNP, pval = exp_fmt$pval.exposure)
  clumped <- tryCatch(
    ieugwasr::ld_clump(clump_input, clump_kb = 1000, clump_r2 = 0.01,
                       clump_p = 1, bfile = ld_ref_dir,
                       plink_bin = plink_exe, pop = "EUR"),
    error = function(e) NULL
  )
  if (is.null(clumped) || nrow(clumped) == 0) return(NULL)
  
  exp_clump <- exp_fmt[exp_fmt$SNP %in% clumped$rsid, ]
  if (nrow(exp_clump) < 1) return(NULL)
  
  out_raw <- gwas_data[SNP %in% exp_clump$SNP, ]
  if (nrow(out_raw) == 0) return(NULL)
  out_fmt <- format_data(
    as.data.frame(out_raw),
    type = "outcome", snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval"
  )
  out_fmt$id.outcome <- "Validation"
  
  mr_data <- harmonise_data(exp_clump, out_fmt)
  if (nrow(mr_data) == 0) return(NULL)
  mr_res <- mr(mr_data, method_list = c("mr_ivw", "mr_wald_ratio"))
  ivw <- mr_res[mr_res$method == "Inverse variance weighted", ]
  if (nrow(ivw) == 0) ivw <- mr_res[mr_res$method == "Wald ratio", ]
  if (nrow(ivw) == 0) return(NULL)
  
  return(data.table(gene = g, b = ivw$b, se = ivw$se, pval = ivw$pval, nsnp = ivw$nsnp))
}

# ---------- 5.2 Batch validation function ----------
run_validation_batch <- function(gene_list, eqtl_data, gwas_data, label, batch_size = 100) {
  if (length(gene_list) == 0) return(NULL)
  
  output_file <- file.path("results", paste0("validation_", label, ".csv"))
  progress_file <- file.path("results", paste0("validation_progress_", label, ".txt"))
  
  if (file.exists(progress_file) && !file.exists(output_file)) {
    cat("⚠️ 进度文件存在但无结果文件，将重新开始。\n")
    file.remove(progress_file)
  }
  
  done_genes <- if (file.exists(progress_file)) readLines(progress_file) else character(0)
  remaining_genes <- setdiff(gene_list, done_genes)
  
  if (length(remaining_genes) == 0) {
    cat("所有验证基因已完成。\n")
    return(if (file.exists(output_file)) fread(output_file) else NULL)
  }
  
  n_batches <- ceiling(length(remaining_genes) / batch_size)
  cat("验证将分为", n_batches, "批，每批最多", batch_size, "个基因\n")
  
  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(remaining_genes))
    batch_genes <- remaining_genes[start_idx:end_idx]
    
    cat(sprintf("\n===== 验证批次 %d/%d =====\n", i, n_batches))
    
    res_list <- list()
    success_genes <- c()
    for (g in batch_genes) {
      res <- tryCatch(validate_gene(g, eqtl_data, gwas_data), error = function(e) NULL)
      if (!is.null(res)) {
        res_list[[length(res_list) + 1]] <- res
        success_genes <- c(success_genes, g)
      }
    }
    
    if (length(res_list) > 0) {
      fwrite(rbindlist(res_list, fill = TRUE), output_file, append = file.exists(output_file))
    }
    
    done_genes <- unique(as.character(c(done_genes, success_genes)))
    writeLines(done_genes, progress_file)
    rm(res_list); gc()
    Sys.sleep(1)
  }
  
  if (file.exists(output_file)) {
    cat("✅ 验证", label, "全部完成。\n")
    return(fread(output_file))
  } else {
    cat("❌ 验证", label, "未产生任何结果。\n")
    return(NULL)
  }
}

# ---------- 5.3 Perform all external validations ----------
# 读取 C 级基因列表（已在第 4 部分生成）
C_fdr_genes <- fread("results/integrated_genes_fdr_main.csv")[grade == "C", gene]
C_nominal   <- integrated[grade == "C", gene]   # 来自名义显著整合表

# AB 级基因验证（使用食道黏膜 eQTL）
cat("\n--- AB 级基因外部验证 ---\n")
val_AB_gwas1 <- run_validation_batch(AB_genes, eqtl_eso, ext_gwas1, "AB_GCST90041888")
val_AB_gwas2 <- run_validation_batch(AB_genes, eqtl_eso, ext_gwas2, "AB_GCST012237")

# FDR‑C 级基因验证（使用唾液腺 eQTL）
cat("\n--- FDR‑C 级基因外部验证 ---\n")
val_Cfdr_gwas1 <- run_validation_batch(C_fdr_genes, eqtl_sal, ext_gwas1, "Cfdr_GCST90041888")
val_Cfdr_gwas2 <- run_validation_batch(C_fdr_genes, eqtl_sal, ext_gwas2, "Cfdr_GCST012237")

# 名义‑C 级基因验证
cat("\n--- 名义‑C 级基因外部验证 ---\n")
val_Cnom_gwas1 <- run_validation_batch(C_nominal, eqtl_sal, ext_gwas1, "Cnom_GCST90041888")
val_Cnom_gwas2 <- run_validation_batch(C_nominal, eqtl_sal, ext_gwas2, "Cnom_GCST012237")

# 排除集基因外部验证 
cat("\n--- 排除集基因外部验证 ---\n")
val_blood_excl_gwas1  <- run_validation_batch(excl_blood,  eqtl_blood,  ext_gwas1, "blood_excluded_GCST90041888")
val_muscle_excl_gwas1 <- run_validation_batch(excl_muscle, eqtl_muscle, ext_gwas1, "muscle_excluded_GCST90041888")

val_blood_excl_gwas2  <- run_validation_batch(excl_blood,  eqtl_blood,  ext_gwas2, "blood_excluded_GCST012237")
val_muscle_excl_gwas2 <- run_validation_batch(excl_muscle, eqtl_muscle, ext_gwas2, "muscle_excluded_GCST012237")

# ---------- 5.4 Perform lever-I/II/IIIgene filtering  ----------
cat("\n===== 筛选优先基因（基于食管黏膜名义显著基因） =====\n")

eso_nominal <- mr_eso[ivw_pval < 0.05, ]

# 合并验证 1 结果（AB 基因用的是食道黏膜，所以直接用 val_AB_gwas1）
candidates <- merge(
  eso_nominal[, .(gene, eso_b = ivw_b, eso_se = ivw_se)],
  val_AB_gwas1[, .(gene, val1_b = b, val1_se = se, val1_p = pval)],
  by = "gene"
)

# 方向一致性
candidates[, direction_concordant := sign(eso_b) == sign(val1_b)]

# CI 重叠判断
candidates[, eso_ci_low := eso_b - 1.96 * eso_se]
candidates[, eso_ci_high := eso_b + 1.96 * eso_se]
candidates[, val1_ci_low := val1_b - 1.96 * val1_se]
candidates[, val1_ci_high := val1_b + 1.96 * val1_se]
candidates[, ci_overlap := (eso_ci_low <= val1_ci_high) & (val1_ci_low <= eso_ci_high)]

# I 级：方向一致 + CI 重叠 + 验证 1 名义显著
priority_I <- candidates[direction_concordant == TRUE & ci_overlap == TRUE & val1_p < 0.05, ]

# II 级：方向一致 + CI 不重叠 + 验证 1 名义显著
priority_II <- candidates[direction_concordant == TRUE & ci_overlap == FALSE & val1_p < 0.05, ]

# III 级：方向不一致 + 验证 1 名义显著
priority_III <- candidates[direction_concordant == FALSE & val1_p < 0.05, ]

# 补充验证 2 信息
priority_I  <- merge(priority_I,  val_AB_gwas2[, .(gene, val2_b = b, val2_p = pval)], by = "gene", all.x = TRUE)
priority_II <- merge(priority_II, val_AB_gwas2[, .(gene, val2_b = b, val2_p = pval)], by = "gene", all.x = TRUE)
priority_III <- merge(priority_III, val_AB_gwas2[, .(gene, val2_b = b, val2_p = pval)], by = "gene", all.x = TRUE)

# 输出各级基因数量
cat("\n=== 各级基因数量 ===\n")
cat("优先I级（方向一致+CI重叠+验证1显著）:", nrow(priority_I), "\n")
cat("优先II级（方向一致+CI不重叠+验证1显著）:", nrow(priority_II), "\n")
cat("优先III级（方向不一致+验证1显著）:", nrow(priority_III), "\n")

# 保存各级基因列表
fwrite(priority_I,  "results/priority_I_genes.csv")
fwrite(priority_II, "results/priority_II_genes.csv")
fwrite(priority_III,"results/priority_III_genes.csv")

# 统计摘要
priority_summary <- data.table(
  Level = c("I", "II", "III"),
  Count = c(nrow(priority_I), nrow(priority_II), nrow(priority_III)),
  Description = c("方向一致+CI重叠+验证1显著", "方向一致+CI不重叠+验证1显著", "方向不一致+验证1显著")
)
fwrite(priority_summary, "results/priority_gene_summary.csv")
print(priority_summary)

# ---------- 5.5 Generate external-validation statistics for table output ----------
cat("\n===== 外部验证统计汇总 =====\n")

summarize_val <- function(val_df, total_genes, set_name) {
  if (is.null(val_df) || nrow(val_df) == 0) {
    return(data.table(GeneSet = set_name, Total = total_genes, Success = 0, NominalSig = 0))
  }
  data.table(
    GeneSet    = set_name,
    Total      = total_genes,
    Success    = nrow(val_df),
    NominalSig = sum(val_df$pval < 0.05, na.rm = TRUE)
  )
}

val_stats <- rbind(
  summarize_val(val_AB_gwas1,      length(AB_genes),    "AB级（食道黏膜）"),
  summarize_val(val_Cfdr_gwas1,    length(C_fdr_genes),  "FDR-C级（唾液腺）"),
  summarize_val(val_Cnom_gwas1,    length(C_nominal),    "名义-C级（唾液腺）"),
  summarize_val(val_blood_excl_gwas1,  length(excl_blood),  "血液排除集"),
  summarize_val(val_muscle_excl_gwas1, length(excl_muscle), "肌肉排除集")
)

fwrite(val_stats, "results/validation_summary.csv")
print(val_stats)
cat("✅ 外部验证全部完成\n")

# ---------- 5.6 Merge PP.H4 ----------
if (file.exists("results/coloc_AB.csv")) {
  coloc_AB <- fread("results/coloc_AB.csv")
  coloc_max <- coloc_AB[, .(PP.H4_max = max(PP.H4, na.rm = TRUE)), by = gene]
  priority_I <- merge(priority_I, coloc_max, by = "gene", all.x = TRUE)
  fwrite(priority_I, "results/priority_I_genes_with_PP.H4.csv")
  cat("✅ 优先I级基因与PP.H4合并完成\n")
}

# ========================== 6. Functional enrichment analysis===================================================

cat("\n===== 6. 功能富集分析 Functional enrichment analysis=====\n")

# ---------- 6.1 Gene-ID conversion function ----------
convert_ids <- function(ens_list) {
  ens_clean <- sub("\\..*", "", ens_list)
  res <- bitr(ens_clean, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Hs.eg.db)
  unique(res$ENTREZID)
}

# ---------- 6.2 Enrichment analysis for grade-A/B genes----------
ab_genes_df <- integrated[grade %in% c("A","B"), .(gene, eso_b)]
ab_genes_df <- na.omit(ab_genes_df)
prot_genes_ab <- ab_genes_df[eso_b < 0, gene]
risk_genes_ab <- ab_genes_df[eso_b > 0, gene]

# 保护性基因
prot_entrez <- convert_ids(prot_genes_ab)
if (length(prot_entrez) > 10) {
  ego_prot <- enrichGO(gene = prot_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                       pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego_prot) && nrow(ego_prot) > 0) {
    fwrite(as.data.frame(ego_prot), "results/GO_BP_Protective.csv")
    pdf("results/figures/GO_Protective_Bubble.pdf", width = 10, height = 8)
    print(dotplot(ego_prot, showCategory = 20, title = "GO BP - Protective Genes"))
    dev.off()
  }
  kegg_prot <- enrichKEGG(gene = prot_entrez, organism = "hsa", pvalueCutoff = 0.05)
  if (!is.null(kegg_prot) && nrow(kegg_prot) > 0) {
    fwrite(as.data.frame(kegg_prot), "results/KEGG_Protective.csv")
    pdf("results/figures/KEGG_Protective_Bubble.pdf", width = 10, height = 6)
    print(dotplot(kegg_prot, showCategory = 15, title = "KEGG - Protective Genes"))
    dev.off()
  }
}

# 风险性基因
risk_entrez <- convert_ids(risk_genes_ab)
if (length(risk_entrez) > 10) {
  ego_risk <- enrichGO(gene = risk_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                       pvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(ego_risk) && nrow(ego_risk) > 0) {
    fwrite(as.data.frame(ego_risk), "results/GO_BP_Risk.csv")
    pdf("results/figures/GO_Risk_Bubble.pdf", width = 10, height = 8)
    print(dotplot(ego_risk, showCategory = 20, title = "GO BP - Risk Genes"))
    dev.off()
  }
  kegg_risk <- enrichKEGG(gene = risk_entrez, organism = "hsa", pvalueCutoff = 0.05)
  if (!is.null(kegg_risk) && nrow(kegg_risk) > 0) {
    fwrite(as.data.frame(kegg_risk), "results/KEGG_Risk.csv")
    pdf("results/figures/KEGG_Risk_Bubble.pdf", width = 10, height = 6)
    print(dotplot(kegg_risk, showCategory = 15, title = "KEGG - Risk Genes"))
    dev.off()
  }
}

# ---------- 6.3 Enrichment analysis for FDR-C grade genes----------
if (file.exists("results/integrated_genes_fdr_main.csv")) {
  integrated_fdr <- fread("results/integrated_genes_fdr_main.csv")
  C_fdr_genes <- integrated_fdr[grade == "C", gene]
  
  if (length(C_fdr_genes) > 10) {
    C_fdr_entrez <- convert_ids(C_fdr_genes)
    if (length(C_fdr_entrez) > 10) {
      # GO 富集
      ego_C <- enrichGO(gene = C_fdr_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                        pvalueCutoff = 0.05, readable = TRUE)
      if (!is.null(ego_C) && nrow(ego_C) > 0) {
        fwrite(as.data.frame(ego_C), "results/GO_BP_C_fdr.csv")
        pdf("results/figures/GO_C_fdr_Bubble.pdf", width = 10, height = 8)
        print(dotplot(ego_C, showCategory = 20, title = "GO BP - FDR-C Genes"))
        dev.off()
      }
      # KEGG 富集
      kegg_C <- enrichKEGG(gene = C_fdr_entrez, organism = "hsa", pvalueCutoff = 0.05)
      if (!is.null(kegg_C) && nrow(kegg_C) > 0) {
        fwrite(as.data.frame(kegg_C), "results/KEGG_C_fdr.csv")
        pdf("results/figures/KEGG_C_fdr_Bubble.pdf", width = 10, height = 6)
        print(dotplot(kegg_C, showCategory = 15, title = "KEGG - FDR-C Genes"))
        dev.off()
      }
    }
  }
} else {
  warning("FDR 整合文件不存在，跳过 FDR-C 富集分析")
}

cat("✅ 富集分析完成\n")

# ========================== Gene Set Enrichment Analysis ==========================
cat("\n===== GSEA 分析 =====\n")
library(clusterProfiler)
library(org.Hs.eg.db)

mr_eso <- fread("results/MR_results_Esophagus_Mucosa.csv")
geneList <- mr_eso$ivw_b
names(geneList) <- sub("\\..*", "", mr_eso$gene)

entrez_map <- bitr(names(geneList), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
idx <- match(names(geneList), entrez_map$ENSEMBL)
keep <- !is.na(idx)
geneList_entrez <- geneList[keep]
names(geneList_entrez) <- entrez_map$ENTREZID[idx[keep]]
geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]
geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

gsea_kegg <- gseKEGG(
  geneList = geneList_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.25,
  verbose = FALSE
)

if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
  fwrite(as.data.frame(gsea_kegg@result), "results/GSEA_KEGG_protective.csv")
  pdf("results/figures/SuppFigure_GSEA_KEGG.pdf", width = 10, height = 6)
  print(dotplot(gsea_kegg, showCategory = 15, title = "GSEA KEGG"))
  dev.off()
}
cat("✅ GSEA 分析完成\n")

# ========================== 7. TCGA oral cancer tissue validation==============================

dir.create("results/TCGA", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# 加载包
pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "edgeR",
          "data.table", "dplyr", "ggplot2", "pheatmap", "ppcor",
          "org.Hs.eg.db", "progress", "BiocManager")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "edgeR", "org.Hs.eg.db"))
      BiocManager::install(pkg, update = FALSE)
    else
      install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

if (!exists("mr_eso")) {
  mr_eso <- fread("results/MR_results_Esophagus_Mucosa.csv")
}

# ---- 7.1. Protective‑gene list ----
protective_genes <- mr_eso[ivw_b < 0 & ivw_pval < 0.05, unique(gene)]
fwrite(data.table(gene = protective_genes), "results/protective_genes_from_MR.csv")
cat("保护性基因数:", length(protective_genes), "\n")

# ---- 7.2. Automatically detect and load TCGA datasets ----
DATA_DIR <- file.path(WORK_DIR, "data")
gdc_cache <- file.path(DATA_DIR, "GDCdata")

query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# 检查 GDCdata 目录下是否已有足够的数据文件
data_exists <- dir.exists(gdc_cache) && 
  length(list.files(gdc_cache, recursive = TRUE, pattern = "star_gene_counts\\.tsv$")) > 500

if (data_exists) {
  cat("✅ 检测到已有完整 TCGA 数据，跳过下载，直接读取。\n")
} else {
  cat("未检测到完整数据，开始下载...\n")
  # 清理可能残留的损坏文件
  if (dir.exists(gdc_cache)) unlink(gdc_cache, recursive = TRUE)
  GDCdownload(query, directory = DATA_DIR, method = "api", files.per.chunk = 10)
}

tcga_data <- tryCatch({
  GDCprepare(query, directory = DATA_DIR, summarizedExperiment = TRUE)
}, error = function(e) {
  GDCprepare(query, directory = DATA_DIR, summarizedExperiment = TRUE, add.gistic2 = FALSE)
})

counts_raw <- assay(tcga_data, "unstranded")
colnames(counts_raw) <- substr(colnames(counts_raw), 1, 15)
clinical <- colData(tcga_data)
cat("✅ 数据读取完成：", nrow(counts_raw), "个基因,", ncol(counts_raw), "个样本\n")

# ---- 7.3. Filter oral cavity subsites ----
site_col <- "tissue_or_organ_of_origin"
oral_keywords <- c("Ventral surface of tongue","Upper Gum","Retromolar area",
                   "Mouth, NOS","Lip, NOS","Hard palate","Gum, NOS",
                   "Floor of mouth","Cheek mucosa")
exclude_keywords <- c("Base of tongue","Tonsil","Oropharynx","Nasopharynx",
                      "Hypopharynx","Larynx","Pharynx","Overlapping")
is_oral <- grepl(paste(oral_keywords, collapse="|"), clinical[[site_col]], ignore.case=TRUE) &
  !grepl(paste(exclude_keywords, collapse="|"), clinical[[site_col]], ignore.case=TRUE)
oral_clinical <- clinical[is_oral, ]
rownames(oral_clinical) <- substr(rownames(oral_clinical), 1, 15)
oral_samples <- colnames(counts_raw)[colnames(counts_raw) %in% rownames(oral_clinical)]
tumor_samples <- oral_samples[grep("-01", oral_samples)]
normal_samples <- oral_samples[grep("-11", oral_samples)]
cat("肿瘤样本:", length(tumor_samples), "正常样本:", length(normal_samples), "\n")

counts_oral <- counts_raw[, c(tumor_samples, normal_samples)]
genes_clean <- sub("\\..*", "", rownames(counts_oral))
counts_oral <- counts_oral[!duplicated(genes_clean), ]
rownames(counts_oral) <- genes_clean[!duplicated(genes_clean)]

# ---- 7.4. LogCPM and protective score ----
dge <- DGEList(counts = counts_oral)
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

common_prot <- intersect(protective_genes, rownames(logcpm))
cat("匹配保护性基因:", length(common_prot), "\n")
if (length(common_prot) < 5) stop("保护性基因数量不足")
prot_score <- colMeans(logcpm[common_prot, , drop = FALSE])

# ---- 7.5. Immune cell signatures combined with ESTIMATE ----
immune_sets <- list(
  CD8_T   = c("CD8A","CD8B","GZMA","GZMB","PRF1"),
  B_cell  = c("CD19","CD79A","CD79B","MS4A1"),
  NK      = c("NKG7","GNLY","KLRD1","KLRF1"),
  M1      = c("NOS2","IL12A","IL12B","TNF","CCL3"),
  M2      = c("MRC1","ARG1","IL10","CD163"),
  DC      = c("CD1C","CD83","CCL17","FCER1A")
)
immune_ens <- lapply(immune_sets, function(gs) {
  ids <- mapIds(org.Hs.eg.db, keys = gs, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  na.omit(ids)
})
immune_scores <- sapply(immune_ens, function(gs) {
  gs_in <- intersect(gs, rownames(logcpm))
  if (length(gs_in) >= 3) colMeans(logcpm[gs_in, , drop = FALSE]) else rep(NA, ncol(logcpm))
})
colnames(immune_scores) <- names(immune_sets)
immune_scores <- immune_scores[, !apply(immune_scores, 2, function(x) all(is.na(x))), drop = FALSE]

est_genes <- c("CD3D","CD3E","CD3G","CD2","CD6","CD7","CD8A","CD8B","GZMA","GZMB",
               "PRF1","IFNG","TNF","LTB","CD19","CD79A","CD79B","MS4A1","TNFRSF17",
               "BLK","PAX5","NKG7","GNLY","KLRD1","KLRF1","KLRB1","NCR1","NCR3",
               "CD14","CD68","CD163","FCGR3A","FCGR3B","CSF1R","CD1C","CLEC10A",
               "CLEC4C","FCER1A","NRP1","HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1",
               "HLA-DPA1","HLA-DPB1")
est_ens <- na.omit(mapIds(org.Hs.eg.db, keys = est_genes, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first"))
common_est <- intersect(est_ens, rownames(logcpm))
immune_score_manual <- colMeans(logcpm[common_est, , drop = FALSE])

# ---- 7.6a. Partial correlation ----
library(ppcor)
partial_res <- data.table()
for (feat in colnames(immune_scores)) {
  rc <- cor.test(prot_score, immune_scores[, feat], method = "spearman")
  pc <- pcor.test(prot_score, immune_scores[, feat], immune_score_manual, method = "spearman")
  partial_res <- rbind(partial_res, data.table(
    Feature = feat, raw_rho = rc$estimate, raw_p = rc$p.value,
    partial_rho = pc$estimate, partial_p = pc$p.value
  ))
}
fwrite(partial_res, "results/TCGA/Partial_Correlation_TCGA.csv")
print(partial_res)

# ---- 7.6b. Permutation test (n=10000) ----
set.seed(2024)
n_perm <- 10000
# 预分配矩阵存储置换后的偏相关系数
perm_rho_matrix <- matrix(NA, nrow = n_perm, ncol = ncol(immune_scores))
colnames(perm_rho_matrix) <- colnames(immune_scores)
for (i in 1:n_perm) {
  # 随机打乱保护性评分
  perm_score <- sample(prot_score)
  for (feat in colnames(immune_scores)) {
    pc <- tryCatch(
      pcor.test(perm_score, immune_scores[, feat], immune_score_manual, method = "spearman"),
      error = function(e) list(estimate = NA)
    )
    perm_rho_matrix[i, feat] <- pc$estimate
  }
}
# 计算观察到的偏相关系数
obs_rho <- sapply(colnames(immune_scores), function(feat) {
  pcor.test(prot_score, immune_scores[, feat], immune_score_manual, method = "spearman")$estimate
})
# 计算经验 P 值（双侧）
empirical_p <- sapply(colnames(immune_scores), function(feat) {
  obs_val <- obs_rho[feat]
  perm_vals <- perm_rho_matrix[, feat]
  (sum(abs(perm_vals) >= abs(obs_val), na.rm = TRUE) + 1) / (n_perm + 1)
})
perm_out <- data.frame(Feature = names(empirical_p), Empirical_P = empirical_p)
fwrite(perm_out, "results/TCGA/Permutation_Test.csv")
cat("✅ 置换检验完成\n")
# ---- 7.7. Differentially expressed genes ----
priority_I <- fread("results/priority_I_genes.csv")
priority_symbols <- mapIds(org.Hs.eg.db, keys = priority_I$gene, 
                           column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
priority_symbols <- priority_symbols[!is.na(priority_symbols)]
priority_ens <- mapIds(org.Hs.eg.db, keys = priority_symbols, 
                       column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
priority_ens <- priority_ens[!is.na(priority_ens)]

counts_sub <- counts_oral[rownames(counts_oral) %in% priority_ens, ]
sample_info <- data.frame(
  row.names = c(tumor_samples, normal_samples),
  group = factor(c(rep("Tumor", length(tumor_samples)), rep("Normal", length(normal_samples))),
                 levels = c("Normal", "Tumor"))
)
dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = sample_info, design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "Tumor", "Normal"))
res <- as.data.frame(res[order(res$padj), ])
diff_table <- data.table(
  Gene = mapIds(org.Hs.eg.db, keys = rownames(res), column = "SYMBOL", keytype = "ENSEMBL"),
  Ensembl = rownames(res), log2FC = res$log2FoldChange, pvalue = res$pvalue, padj = res$padj
)
fwrite(diff_table, "results/TCGA/Table_DiffExpression_11Genes_DESeq2.csv")
cat("✅ 11基因差异表达分析完成\n")

# ---- 7.8. Box plot ----
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
expr_norm <- assay(vsd)
pdf("results/figures/Figure_Boxplot_11Genes.pdf", width = 14, height = 8)
par(mfrow = c(3, 4))
for (g in rownames(expr_norm)) {
  gene_sym <- mapIds(org.Hs.eg.db, keys = g, column = "SYMBOL", keytype = "ENSEMBL")
  boxplot(expr_norm[g, ] ~ sample_info$group, main = gene_sym, ylab = "VST expression", xlab = "")
}
dev.off()
cat("✅ 箱线图已保存\n")

# ---- 7.9. Heatmap ----
if (exists("prot_score") && exists("immune_scores")) {
  cor_mat <- cor(cbind(prot_score, immune_scores), method = "spearman", use = "pairwise.complete.obs")
  pheatmap(cor_mat, display_numbers = TRUE, number_format = "%.2f",
           color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
           main = "Protective Score vs Immune Features",
           filename = "results/figures/Figure3_Heatmap.pdf", width = 8, height = 7)
  cat("✅ 热图已保存\n")
}

cat("\n========== TCGA 组织验证全部完成！ ==========\n")


# ---- 7.10 免疫谱系标志基因重叠分析（使用全部 229 个保护性基因） ----
cat("\n===== 免疫谱系重叠分析 =====\n")

# 严格按照论文补充表 S1 的基因列表
immune_lineage_genes <- list(
  B_cells         = c("CD19","CD79A","CD79B","MS4A1","PAX5","BLK"),
  CD4_T_cells     = c("CD4"),
  CD8_T_cells     = c("CD8A","CD8B"),
  NK_cells        = c("NKG7","GNLY","KLRD1","KLRF1","NCR1","NCR3"),
  Monocytes       = c("CD14","FCGR3A","CSF1R"),
  Dendritic_cells = c("CLEC4C","NR4A1","CLEC10A","FCER1A"),
  Macrophages     = c("CD68","CD163","MRC1")
)

immune_lineage_ens <- lapply(immune_lineage_genes, function(gs) {
  ids <- mapIds(org.Hs.eg.db, keys = gs, keytype = "SYMBOL", 
                column = "ENSEMBL", multiVals = "first")
  na.omit(ids)
})

# 注意：这里使用全部 229 个保护性基因，而非 TCGA 匹配后的子集
overlap_results <- data.frame(
  CellType       = names(immune_lineage_ens),
  LineageGenes   = sapply(immune_lineage_ens, length),
  Overlap        = sapply(immune_lineage_ens, function(x) length(intersect(x, protective_genes))),
  OverlapPercent = sapply(immune_lineage_ens, function(x) 
    round(length(intersect(x, protective_genes)) / length(x) * 100, 1))
)

fwrite(overlap_results, "results/TCGA/Table_CellType_Specificity.csv", row.names = FALSE)
print(overlap_results)
cat("✅ 免疫谱系重叠分析完成\n")

# ---- 7.11 风险评分对照分析 + 风险评分置换检验 ----
cat("\n===== 风险评分对照分析 =====\n")

risk_genes <- mr_eso[ivw_b > 0 & ivw_pval < 0.05, unique(gene)]
common_risk <- intersect(risk_genes, rownames(logcpm))

if (length(common_risk) >= 5) {
  # 计算风险基因评分
  risk_score <- colMeans(logcpm[common_risk, , drop = FALSE])
  
  # 风险评分与免疫特征的相关性
  risk_cor_results <- data.frame()
  for (feat in colnames(immune_scores)) {
    rc <- cor.test(risk_score, immune_scores[, feat], method = "spearman", use = "complete.obs")
    risk_cor_results <- rbind(risk_cor_results, data.frame(
      Feature = feat, Risk_rho = rc$estimate, Risk_p = rc$p.value
    ))
  }
  fwrite(risk_cor_results, "results/TCGA/Risk_Score_Correlation.csv", row.names = FALSE)
  
  # 保护性评分与风险评分对比
  prot_cor <- data.frame()
  for (feat in colnames(immune_scores)) {
    rc <- cor.test(prot_score, immune_scores[, feat], method = "spearman", use = "complete.obs")
    prot_cor <- rbind(prot_cor, data.frame(
      Feature = feat, Protective_rho = rc$estimate, Protective_p = rc$p.value
    ))
  }
  comparison <- merge(prot_cor, risk_cor_results, by = "Feature")
  fwrite(comparison, "results/TCGA/Protective_vs_Risk_Comparison.csv", row.names = FALSE)
  
  # 风险评分置换检验
  set.seed(2024)
  n_perm <- 10000
  perm_risk_cor <- matrix(NA, nrow = n_perm, ncol = ncol(immune_scores))
  colnames(perm_risk_cor) <- colnames(immune_scores)
  for (i in 1:n_perm) {
    rand_genes <- sample(rownames(logcpm), length(common_risk))
    rand_score <- colMeans(logcpm[rand_genes, , drop = FALSE])
    for (feat in colnames(immune_scores)) {
      perm_risk_cor[i, feat] <- cor(rand_score, immune_scores[, feat], 
                                    method = "spearman", use = "complete.obs")
    }
  }
  risk_empirical_p <- sapply(colnames(immune_scores), function(feat) {
    obs <- risk_cor_results$Risk_rho[risk_cor_results$Feature == feat]
    perm_vals <- perm_risk_cor[, feat]
    if (obs >= 0) {
      (sum(perm_vals >= obs, na.rm = TRUE) + 1) / (n_perm + 1)
    } else {
      (sum(perm_vals <= obs, na.rm = TRUE) + 1) / (n_perm + 1)
    }
  })
  fwrite(data.frame(Feature = names(risk_empirical_p), Empirical_P = risk_empirical_p),
         "results/TCGA/Risk_Score_Permutation.csv", row.names = FALSE)
  
  cat("✅ 风险评分对照及置换检验完成\n")
} else {
  cat("⚠️ 风险基因数量不足（<5），跳过风险评分对照分析\n")
}

# ---- 7.12 留一法敏感性分析（CD8/NK） ----
cat("\n===== 留一法敏感性分析 =====\n")

cd8_ens <- na.omit(mapIds(org.Hs.eg.db, keys = c("CD8A","CD8B"), 
                          column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first"))
nk_ens  <- na.omit(mapIds(org.Hs.eg.db, keys = c("NKG7","GNLY","KLRD1","KLRF1","NCR1","NCR3"), 
                          column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first"))

sensitivity_res <- data.table()
for (feat in c("CD8_T", "NK")) {
  for (adj in c("Full","No_CD8","No_NK","No_Both")) {
    if (adj == "Full") {
      score_adj <- immune_score_manual
    } else if (adj == "No_CD8") {
      score_adj <- colMeans(logcpm[setdiff(common_est, cd8_ens), , drop = FALSE])
    } else if (adj == "No_NK") {
      score_adj <- colMeans(logcpm[setdiff(common_est, nk_ens), , drop = FALSE])
    } else {
      score_adj <- colMeans(logcpm[setdiff(common_est, union(cd8_ens, nk_ens)), , drop = FALSE])
    }
    pc <- pcor.test(prot_score, immune_scores[, feat], score_adj, method = "spearman")
    sensitivity_res <- rbind(sensitivity_res, data.table(
      Feature = feat, Adjustment = adj, 
      partial_rho = pc$estimate, partial_p = pc$p.value
    ))
  }
}
fwrite(sensitivity_res, "results/TCGA/LeaveOneOut_Sensitivity.csv", row.names = FALSE)
print(sensitivity_res)
cat("✅ 留一法敏感性分析完成\n")

# ---- 7.13 Bootstrap 置信区间（Spearman ρ） ----
cat("\n===== Bootstrap 置信区间计算 =====\n")

set.seed(123)
bootstrap_results <- data.frame(Feature = character(), rho = numeric(),
                                CI_lower = numeric(), CI_upper = numeric())
for (feat in colnames(immune_scores)) {
  data_boot <- data.frame(prot = prot_score, imm = immune_scores[, feat])
  data_boot <- na.omit(data_boot)
  boot_func <- function(data, indices) cor(data[indices,1], data[indices,2], method = "spearman")
  boot_res <- boot(data_boot, boot_func, R = 1000)
  ci <- boot.ci(boot_res, type = "perc")
  bootstrap_results <- rbind(bootstrap_results, data.frame(
    Feature = feat,
    rho = cor(data_boot$prot, data_boot$imm, method = "spearman"),
    CI_lower = if (!is.null(ci)) ci$percent[4] else NA,
    CI_upper = if (!is.null(ci)) ci$percent[5] else NA
  ))
}
fwrite(bootstrap_results, "results/TCGA/Table_Spearman_CI.csv", row.names = FALSE)
print(bootstrap_results)
cat("✅ Bootstrap 置信区间计算完成\n")
# ========================== 8. Supplementary figures==============================

library(ggplot2)
library(data.table)
library(ggrepel)

cat("\n===== 8. Supplementary figures =====\n")
# Q-Q 图
mr_eso <- fread("results/MR_results_Esophagus_Mucosa.csv")
observed_p <- sort(mr_eso$ivw_pval)
expected_p <- ppoints(length(observed_p))
qq_df <- data.frame(Expected = -log10(expected_p), Observed = -log10(observed_p))
ggplot(qq_df, aes(x = Expected, y = Observed)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "-log10(Expected P)", y = "-log10(Observed P)", title = "Q-Q Plot") + theme_bw()
ggsave("results/figures/SuppFigure_QQ.pdf", width = 6, height = 6)

#Figure1A_Volcano
integrated <- fread("results/integrated_genes_nominal_exploratory.csv")
volcano <- merge(mr_eso, integrated[, .(gene, grade)], by = "gene", all.x = TRUE)
volcano[, log10p := -log10(ivw_pval)]
volcano[, col_group := ifelse(grade %in% c("A","B"), grade, "Other")]
priority_I <- fread("results/priority_I_genes.csv")
volcano[, label := ifelse(gene %in% priority_I$gene,
                          mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"), "")]
ggplot(volcano, aes(x = ivw_b, y = log10p, color = col_group)) + geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("A" = "#E41A1C", "B" = "#377EB8", "Other" = "grey70")) +
  labs(x = "MR Effect (β)", y = "-log10(P)", title = "Esophageal Mucosa MR Results") + theme_minimal()
ggsave("results/figures/Figure1A_Volcano.pdf", width = 7, height = 5)
cat("✅ 补充图表完成\n")

#Figure 1B
mr_sal <- fread("results/MR_results_Minor_Salivary_Gland.csv")
sal_plot <- merge(mr_sal, integrated[, .(gene, grade)], by = "gene", all.x = TRUE)
sal_plot[, log10p := -log10(ivw_pval)]
sal_plot[, col_group := ifelse(grade == "C", "FDR-C", "Other")]

ggplot(sal_plot, aes(x = ivw_b, y = log10p, color = col_group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("FDR-C" = "#FF7F00", "Other" = "grey70")) +
  labs(x = "MR Effect (β)", y = "-log10(P)", 
       title = "Minor Salivary Gland (FDR-C)") +
  theme_minimal()
ggsave("results/figures/Figure1B_Salivary.pdf", width = 7, height = 5)
cat("✅ Figure1B 唾液腺散点图已保存\n")

#森林图 (优先 I 级基因 OR 95% CI)
priority_I <- fread("results/priority_I_genes.csv")

# 获取基因Symbol
priority_I[, Symbol := mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL", 
                              keytype = "ENSEMBL", multiVals = "first")]
# 计算 OR 和 95% CI
priority_I[, `:=`(OR = exp(eso_b),
                  OR_low = exp(eso_b - 1.96 * eso_se),
                  OR_high = exp(eso_b + 1.96 * eso_se))]

# 按 OR 排序
priority_I <- priority_I[order(OR)]

# 绘制森林图
ggplot(priority_I, aes(x = OR, y = reorder(Symbol, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  scale_x_log10() +
  labs(x = "Odds Ratio (95% CI)", y = "", 
       title = "Priority I Genes - Forest Plot") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("results/figures/Figure2_Forest.pdf", width = 8, height = 6)
cat("✅ 森林图已保存\n")

# 保护性基因 GO 富集气泡图 
if (file.exists("results/GO_BP_Protective.csv")) {
  go_prot <- fread("results/GO_BP_Protective.csv")
  go_prot[, GeneRatio := sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))]
  go_top20 <- head(go_prot[order(pvalue), ], 20)
  
  ggplot(go_top20, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                       size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    scale_size_continuous(name = "Gene count") +
    labs(x = "Gene Ratio", y = "", title = "GO BP - Protective Genes") +
    theme_bw()
  ggsave("results/figures/GO_Protective_Bubble.pdf", width = 10, height = 8)
  cat("✅ 保护性基因 GO 气泡图已保存\n")
}

#保护性基因 KEGG 富集气泡图
if (file.exists("results/KEGG_Protective.csv")) {
  kegg_prot <- fread("results/KEGG_Protective.csv")
  kegg_prot[, GeneRatio := sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))]
  kegg_top15 <- head(kegg_prot[order(pvalue), ], 15)
  
  ggplot(kegg_top15, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                         size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    scale_size_continuous(name = "Gene count") +
    labs(x = "Gene Ratio", y = "", title = "KEGG - Protective Genes") +
    theme_bw()
  ggsave("results/figures/KEGG_Protective_Bubble.pdf", width = 10, height = 6)
  cat("✅ 保护性基因 KEGG 气泡图已保存\n")
}

# 风险性基因 GO 富集气泡图
if (file.exists("results/GO_BP_Risk.csv")) {
  go_risk <- fread("results/GO_BP_Risk.csv")
  go_risk[, GeneRatio := sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))]
  go_risk_top20 <- head(go_risk[order(pvalue), ], 20)
  
  ggplot(go_risk_top20, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                            size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    scale_size_continuous(name = "Gene count") +
    labs(x = "Gene Ratio", y = "", title = "GO BP - Risk Genes") +
    theme_bw()
  ggsave("results/figures/GO_Risk_Bubble.pdf", width = 10, height = 8)
  cat("✅ 风险性基因 GO 气泡图已保存\n")
}

# 风险性基因 KEGG 富集气泡图 
if (file.exists("results/KEGG_Risk.csv")) {
  kegg_risk <- fread("results/KEGG_Risk.csv")
  kegg_risk[, GeneRatio := sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))]
  kegg_risk_top15 <- head(kegg_risk[order(pvalue), ], 15)
  
  ggplot(kegg_risk_top15, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                              size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    scale_size_continuous(name = "Gene count") +
    labs(x = "Gene Ratio", y = "", title = "KEGG - Risk Genes") +
    theme_bw()
  ggsave("results/figures/KEGG_Risk_Bubble.pdf", width = 10, height = 6)
  cat("✅ 风险性基因 KEGG 气泡图已保存\n")
}

# FDR‑C 级基因 GO 富集气泡图 
if (file.exists("results/GO_BP_C_fdr.csv")) {
  go_C <- fread("results/GO_BP_C_fdr.csv")
  go_C[, GeneRatio := sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))]
  go_C_top20 <- head(go_C[order(pvalue), ], 20)
  
  ggplot(go_C_top20, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                         size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    scale_size_continuous(name = "Gene count") +
    labs(x = "Gene Ratio", y = "", title = "GO BP - FDR-C Genes") +
    theme_bw()
  ggsave("results/figures/GO_C_fdr_Bubble.pdf", width = 10, height = 8)
  cat("✅ FDR-C 基因 GO 气泡图已保存\n")
}

#FDR‑C 级基因 KEGG 富集气泡图 
if (file.exists("results/KEGG_C_fdr.csv")) {
  kegg_C <- fread("results/KEGG_C_fdr.csv")
  kegg_C[, GeneRatio := sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))]
  kegg_C_top15 <- head(kegg_C[order(pvalue), ], 15)
  
  ggplot(kegg_C_top15, aes(x = GeneRatio, y = reorder(Description, GeneRatio),
                           size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue", name = "FDR") +
    scale_size_continuous(name = "Gene count") +
    labs(x = "Gene Ratio", y = "", title = "KEGG - FDR-C Genes") +
    theme_bw()
  ggsave("results/figures/KEGG_C_fdr_Bubble.pdf", width = 10, height = 6)
  cat("✅ FDR-C 基因 KEGG 气泡图已保存\n")
}

cat("\n========== 所有补充图表生成完毕 ==========\n")
# ==========================Full‑dataset checking checklist for manuscript==============================
# 逐项读取结果文件，打印关键统计量，方便对照论文检查

cat("\n========================================\n")
cat("     论文全部数据核对清单\n")
cat("========================================\n")

# ========== Part 1: Core MR findings ==========
cat("\n========== 1. 核心 MR 发现 ==========\n")

cat("\n▶ 食道黏膜 Tier 1:\n")
mr_eso <- fread("results/MR_results_Esophagus_Mucosa.csv")
eso_nominal <- mr_eso[ivw_pval < 0.05, ]
eso_fdr <- p.adjust(mr_eso$ivw_pval, method = "fdr")
cat("  - 总基因数:", nrow(mr_eso), "\n")
cat("  - 名义显著基因数 (P<0.05):", nrow(eso_nominal), "\n")
cat("  - FDR 最小调整 P 值:", format(min(eso_fdr, na.rm = TRUE), digits = 3), "\n")
cat("  - 通过 FDR<0.05 的基因数:", sum(eso_fdr < 0.05, na.rm = TRUE), "\n")

cat("\n▶ 全血 Tier 1:\n")
mr_blood <- fread("results/MR_results_Whole_Blood.csv")
blood_fdr <- p.adjust(mr_blood$ivw_pval, method = "fdr")
cat("  - 总基因数:", nrow(mr_blood), "\n")
cat("  - 名义显著基因数 (P<0.05):", sum(mr_blood$ivw_pval < 0.05), "\n")
cat("  - 通过 FDR<0.05 的基因数:", sum(blood_fdr < 0.05, na.rm = TRUE), "\n")

cat("\n▶ 骨骼肌 Tier 1:\n")
mr_muscle <- fread("results/MR_results_Muscle_Skeletal.csv")
muscle_fdr <- p.adjust(mr_muscle$ivw_pval, method = "fdr")
cat("  - 总基因数:", nrow(mr_muscle), "\n")
cat("  - 名义显著基因数 (P<0.05):", sum(mr_muscle$ivw_pval < 0.05), "\n")
cat("  - 通过 FDR<0.05 的基因数:", sum(muscle_fdr < 0.05, na.rm = TRUE), "\n")

cat("\n▶ 唾液腺 Tier 2 (FDR-C 级):\n")
mr_sal <- fread("results/MR_results_Minor_Salivary_Gland.csv")
sal_fdr <- p.adjust(mr_sal$ivw_pval, method = "fdr")
cat("  - 总基因数:", nrow(mr_sal), "\n")
cat("  - 名义显著基因数 (P<0.05):", sum(mr_sal$ivw_pval < 0.05), "\n")
cat("  - 通过 FDR<0.05 的基因数:", sum(sal_fdr < 0.05, na.rm = TRUE), "\n")
cat("  - FDR-C 级基因数 (论文报告): 180\n")

# ========== Part 2: Evidence integration ==========
cat("\n========== 2. 证据整合 ==========\n")
integrated <- fread("results/integrated_genes_nominal_exploratory.csv")
cat("\n▶ A/B/C 分级统计:\n")
print(table(integrated$grade))

# ========== Part 3: Prioritized gene grading ==========
cat("\n========== 3. 优先基因分级 ==========\n")
if (file.exists("results/priority_gene_summary.csv")) {
  priority_summary <- fread("results/priority_gene_summary.csv")
  cat("\n▶ I/II/III 级基因数量:\n")
  print(priority_summary)
}

cat("\n▶ 优先 I 级基因完整列表 (11 个):\n")
if (file.exists("results/priority_I_genes.csv")) {
  priority_I <- fread("results/priority_I_genes.csv")
  cat("  - 总数:", nrow(priority_I), "\n")
  cat("  - 基因列表:\n")
  for (i in 1:nrow(priority_I)) {
    gene_id <- priority_I$gene[i]
    # 尝试获取基因Symbol
    gene_sym <- tryCatch(
      mapIds(org.Hs.eg.db, keys = gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
      error = function(e) gene_id
    )
    val2_p <- if (is.na(priority_I$val2_p[i])) "NA" else format(priority_I$val2_p[i], digits = 3)
    cat(sprintf("    %2d. %s (%s) | 验证1 P=%s | 验证2 P=%s\n", 
                i, gene_sym, gene_id, 
                format(priority_I$val1_p[i], digits = 3), val2_p))
  }
  
  cat("\n  - 验证二显著的基因 (P<0.05):\n")
  val2_sig <- priority_I[!is.na(val2_p) & val2_p < 0.05, ]
  for (i in 1:nrow(val2_sig)) {
    gene_sym <- tryCatch(
      mapIds(org.Hs.eg.db, keys = val2_sig$gene[i], column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
      error = function(e) val2_sig$gene[i]
    )
    cat(sprintf("    %s (P=%s)\n", gene_sym, format(val2_sig$val2_p[i], digits = 3)))
  }
}

# ========== Part 4: Colocalization analysis ==========
cat("\n========== 4. 共定位分析 ==========\n")
if (file.exists("results/Table3_Summary.csv")) {
  table3 <- fread("results/Table3_Summary.csv")
  cat("\n▶ 共定位统计 (Table 3):\n")
  print(table3)
}

if (file.exists("results/coloc_AB.csv")) {
  coloc_AB <- fread("results/coloc_AB.csv")
  cat("\n▶ AB 级基因共定位:\n")
  cat("  - 成功分析数:", nrow(coloc_AB), "\n")
  cat("  - PP.H4 ≥ 0.5 的基因数:", sum(coloc_AB$PP.H4 >= 0.5, na.rm = TRUE), "\n")
  cat("  - PP.H4 > 0.8 的基因数:", sum(coloc_AB$PP.H4 > 0.8, na.rm = TRUE), "\n")
}

# ========== Part 5: External validation ==========
cat("\n========== 5. 外部验证 ==========\n")
if (file.exists("results/validation_summary.csv")) {
  val_stats <- fread("results/validation_summary.csv")
  cat("\n▶ 外部验证统计:\n")
  print(val_stats)
}

cat("\n▶ AB 级基因验证详情:\n")
if (file.exists("results/validation_AB_GCST90041888.csv")) {
  val1 <- fread("results/validation_AB_GCST90041888.csv")
  cat("  - 验证1 (GCST90041888): 成功", nrow(val1), "个, 显著", sum(val1$pval < 0.05), "个\n")
}
if (file.exists("results/validation_AB_GCST012237.csv")) {
  val2 <- fread("results/validation_AB_GCST012237.csv")
  cat("  - 验证2 (GCST012237): 成功", nrow(val2), "个, 显著", sum(val2$pval < 0.05), "个\n")
}

# ========== Part 6: TCGA tissue validation for oral cancer ==========
cat("\n========== 6. TCGA 组织验证 ==========\n")
if (file.exists("results/TCGA/Partial_Correlation_TCGA.csv")) {
  partial_res <- fread("results/TCGA/Partial_Correlation_TCGA.csv")
  cat("\n▶ 偏相关分析 (Table 4):\n")
  cat(sprintf("%-12s %8s %10s %8s %10s\n", "特征", "原始ρ", "原始P", "偏ρ", "偏P"))
  for (i in 1:nrow(partial_res)) {
    cat(sprintf("%-12s %8.3f %10s %8.3f %10s\n",
                partial_res$Feature[i],
                partial_res$raw_rho[i],
                format(partial_res$raw_p[i], digits = 3),
                partial_res$partial_rho[i],
                format(partial_res$partial_p[i], digits = 3)))
  }
}

if (file.exists("results/TCGA/Permutation_Test.csv")) {
  perm <- fread("results/TCGA/Permutation_Test.csv")
  cat("\n▶ 置换检验 P 值:\n")
  for (i in 1:nrow(perm)) {
    cat(sprintf("  %-12s: %s\n", perm$Feature[i], format(perm$Empirical_P[i], digits = 3)))
  }
}

if (file.exists("results/TCGA/Table_DiffExpression_11Genes_DESeq2.csv")) {
  diff <- fread("results/TCGA/Table_DiffExpression_11Genes_DESeq2.csv")
  cat("\n▶ 11 基因差异表达 (FDR 显著):\n")
  diff_sig <- diff[padj < 0.05, ]
  if (nrow(diff_sig) > 0) {
    for (i in 1:nrow(diff_sig)) {
      cat(sprintf("  %s: log2FC=%.2f, FDR=%s\n", 
                  diff_sig$Gene[i], diff_sig$log2FC[i], format(diff_sig$padj[i], digits = 3)))
    }
  } else {
    cat("  无基因通过 FDR<0.05 (正常, n=4 正常样本)\n")
  }
}

# ========== Part 7: Functional enrichment analysis ==========
cat("\n========== 7. Functional enrichment analysis ==========\n")
kegg_files <- c("results/KEGG_Protective.csv", "results/KEGG_C_fdr.csv")
for (f in kegg_files) {
  if (file.exists(f)) {
    kegg <- fread(f)
    cat(sprintf("\n▶ %s (前 5 条):\n", basename(f)))
    top5 <- head(kegg[order(pvalue), ], 5)
    for (i in 1:nrow(top5)) {
      cat(sprintf("  %s (P=%s, FDR=%s)\n", 
                  top5$Description[i], 
                  format(top5$pvalue[i], digits = 3),
                  format(top5$p.adjust[i], digits = 3)))
    }
  }
}

# ========== Part 8: Data‑completeness inspection ==========
cat("\n========== 8. 数据完整性自检 ==========\n")
required_files <- c(
  "results/MR_results_Esophagus_Mucosa.csv",
  "results/MR_results_Whole_Blood.csv",
  "results/MR_results_Muscle_Skeletal.csv",
  "results/MR_results_Minor_Salivary_Gland.csv",
  "results/integrated_genes_nominal_exploratory.csv",
  "results/integrated_genes_fdr_main.csv",
  "results/priority_I_genes.csv",
  "results/priority_gene_summary.csv",
  "results/Table3_Summary.csv",
  "results/validation_summary.csv",
  "results/TCGA/Partial_Correlation_TCGA.csv",
  "results/TCGA/Permutation_Test.csv",
  "results/TCGA/Table_DiffExpression_11Genes_DESeq2.csv",
  "results/KEGG_Protective.csv",
  "results/KEGG_C_fdr.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  cat("❌ 以下文件缺失:\n")
  for (f in missing_files) cat("  -", f, "\n")
} else {
  cat("✅ 所有必需文件均已生成。\n")
}

cat("\n========================================\n")
cat("     核对完成 — 请逐项对照论文\n")
cat("========================================\n")


# ==========================Supplementary output: key statistical data underlying figures=====================
library(data.table)
library(ggplot2)

# ---- 1. Lambda (λ) value of Q‑Q plot (genomic inflation factor) ----
cat("\n========== QQ 图统计数据 ==========\n")
mr_eso <- fread("results/MR_results_Esophagus_Mucosa.csv")
chisq <- qchisq(mr_eso$ivw_pval, df = 1, lower.tail = FALSE)
lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
cat("QQ 图 λ 值 (基因组膨胀因子):", round(lambda, 4), "\n")
cat("(λ 接近 1 表示无系统性偏差)\n")

# ---- 2. Count of significant genes in volcano plot ----
cat("\n========== 火山图统计数据 ==========\n")
integrated <- fread("results/integrated_genes_nominal_exploratory.csv")
volcano_data <- merge(mr_eso, integrated[, .(gene, grade)], by = "gene", all.x = TRUE)
volcano_data[, significant := ivw_pval < 0.05]

cat("总基因数:", nrow(volcano_data), "\n")
cat("名义显著基因数 (P<0.05):", sum(volcano_data$significant), "\n")
cat("A 级基因数:", sum(volcano_data$grade == "A", na.rm = TRUE), "\n")
cat("B 级基因数:", sum(volcano_data$grade == "B", na.rm = TRUE), "\n")

# 标记的优先 I 级基因
priority_I <- fread("results/priority_I_genes.csv")
volcano_data[, is_priority := gene %in% priority_I$gene]
cat("图上标注的优先 I 级基因数:", sum(volcano_data$is_priority), "\n")

# ---- 3. Range of effect sizes in forest plot ----
cat("\n========== 森林图统计数据 ==========\n")
priority_I <- fread("results/priority_I_genes.csv")
cat("优先 I 级基因 OR 范围:", round(min(exp(priority_I$eso_b), na.rm = TRUE), 2), 
    "-", round(max(exp(priority_I$eso_b), na.rm = TRUE), 2), "\n")
cat("显著性最强的基因:", 
    priority_I[which.min(val1_p), 
               paste0(mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
                      " (P=", format(val1_p, digits = 3), ")")], "\n")

# ---- 4. Correlation‑coefficient matrix for heatmap ----
cat("\n========== 热图统计数据 ==========\n")
if (exists("prot_score") && exists("immune_scores")) {
  cor_mat <- cor(cbind(prot_score, immune_scores), method = "spearman", use = "pairwise.complete.obs")
  cat("保护性评分与免疫特征相关系数矩阵:\n")
  print(round(cor_mat, 3))
} else {
  cat("prot_score 或 immune_scores 变量不存在，请先运行 TCGA 分析。\n")
}

# ---- 5. Sample‑size confirmation for box plots ----
cat("\n========== 箱线图样本数 ==========\n")
cat("肿瘤样本数:", length(tumor_samples), "\n")
cat("正常样本数:", length(normal_samples), "\n")

# ---- 6. Top‑ranked pathways from enrichment bubble plot ----
cat("\n========== 富集分析关键通路 ==========\n")
kegg_prot <- fread("results/KEGG_Protective.csv")
kegg_top5 <- head(kegg_prot[order(pvalue), ], 5)
cat("保护性基因 KEGG 富集 Top 5:\n")
for (i in 1:nrow(kegg_top5)) {
  cat(sprintf("  %d. %s (P=%s, FDR=%s, %d genes)\n", 
              i, kegg_top5$Description[i], 
              format(kegg_top5$pvalue[i], digits = 3),
              format(kegg_top5$p.adjust[i], digits = 3),
              kegg_top5$Count[i]))
}

cat("\n========== 所有统计数据输出完成 ==========\n")

sink("sessionInfo.txt")
sessionInfo()
sink()

# ========================== 9. Final data validation ----------------------------
cat("\n========== 9. 最终数据核对清单 ==========\n")
cat("食道黏膜总基因数:", nrow(mr_eso), "\n")
cat("名义显著基因数 (P<0.05):", eso_nominal[,.N], "\n")
cat("优先 I 级基因数:", nrow(priority_I), "\n")
cat("优先 II 级基因数:", nrow(priority_II), "\n")
cat("优先 III 级基因数:", nrow(priority_III), "\n")
cat("共定位 AB 级支持数:", sum(coloc_AB$PP.H4 >= 0.5, na.rm = TRUE), "\n")
cat("外部验证 AB 级显著数:", sum(val_AB_gwas1$pval < 0.05, na.rm = TRUE), "\n")
cat("TCGA 偏相关 B_cell 偏ρ:", partial_res[Feature == "B_cell", round(partial_rho, 3)], "P:", format(partial_res[Feature == "B_cell", partial_p], digits = 3), "\n")
cat("TCGA 偏相关 CD8_T 偏ρ:", partial_res[Feature == "CD8_T", round(partial_rho, 3)], "P:", format(partial_res[Feature == "CD8_T", partial_p], digits = 3), "\n")

