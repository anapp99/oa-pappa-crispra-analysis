# =====================================================================
# config.R
# Shared settings for all analysis scripts.
# This is the ONLY file containing paths - edit here, nowhere else.
# All other paths are relative to the repository root.
# =====================================================================

# ---------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------

data_dir    <- "data"
results_dir <- "results"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------
# Salmon quantification directory (script 01 only)
#
# Set this to the directory containing one sub-directory per sample,
# each holding a quant.sf file. This is the only machine-specific path
# in the repository and is not required for scripts 02-04, which start
# from the differential expression tables in data/.
# ---------------------------------------------------------------------

salmon_dir <- Sys.getenv("SALMON_DIR", unset = file.path(data_dir, "salmon"))
meta_file  <- file.path(data_dir, "PAPPA_RNAseq_Meta.csv")

# ---------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------

file_crispra_de <- file.path(data_dir, "DE_PAPPA_CRISPRa_vs_SCR.csv")
file_dmc_de     <- file.path(data_dir, "DE_DMC_vs_PLC_paired.csv")
file_hallmark   <- file.path(data_dir, "hallmark_sets.tsv")

# ---------------------------------------------------------------------
# Analysis thresholds
# ---------------------------------------------------------------------

padj_cutoff <- 0.05

# Figure 6A (volcano) and the headline DEG count
lfc_cutoff_volcano <- 1.5

# DESeq2 pre-filter: mean raw count across all samples
min_mean_count <- 5

# ---------------------------------------------------------------------
# GSEA settings
#
# rank_metric: "log2FoldChange" produced the results presented in the
# paper. "stat" (DESeq2 Wald statistic) is run alongside as a
# sensitivity analysis; the leading pathways are unchanged.
# ---------------------------------------------------------------------

gsea_rank_metrics <- c("log2FoldChange", "stat")
gsea_primary_metric <- "log2FoldChange"

gsea_seeds    <- c(123, 1, 42)   # 123 is primary; 1 and 42 are stability checks
gsea_min_size <- 10
gsea_max_size <- 500

# ---------------------------------------------------------------------
# Helper: write a CSV without ever overwriting an existing file
# ---------------------------------------------------------------------

write_csv_safe <- function(x, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop(
      "Refusing to overwrite an existing file:\n  ", path,
      "\nDelete it first, or call with overwrite = TRUE.",
      call. = FALSE
    )
  }
  utils::write.csv(x, path, row.names = FALSE)
  message("Wrote: ", path)
  invisible(path)
}

# ---------------------------------------------------------------------
# Helper: record the session for reproducibility
# ---------------------------------------------------------------------

save_session_info <- function(tag) {
  path <- file.path(results_dir, paste0("sessionInfo_", tag, ".txt"))
  writeLines(capture.output(sessionInfo()), path)
  message("Wrote: ", path)
  invisible(path)
}
