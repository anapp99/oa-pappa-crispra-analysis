# =====================================================================
# 02_gsea_hallmark.R
#
# Pre-ranked Hallmark GSEA for both comparisons:
#   - CRISPRa vs SCR          (data/DE_PAPPA_CRISPRa_vs_SCR.csv)
#   - DMC vs PLC cartilage    (data/DE_DMC_vs_PLC_paired.csv)
#
# Genes are ranked by log2FoldChange (the metric used for the results
# presented in the paper) and, as a sensitivity analysis, by the DESeq2
# Wald statistic. Each is run under three random seeds to confirm the
# results are stable.
#
# Hallmark gene sets are read from a frozen local file
# (data/hallmark_sets.tsv) rather than queried live, so the analysis
# does not change when MSigDB is updated. See R/make_hallmark_sets.R.
#
# Output: results/GSEA_Hallmark_<comparison>_rank<metric>_allseeds.csv
#         results/GSEA_Hallmark_<comparison>_rank<metric>.csv  (seed 123)
#         results/GSEA_rank_metric_comparison.csv
#
# Run from the repository root:  Rscript R/02_gsea_hallmark.R
# =====================================================================

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)   # GSEA() wraps fgsea
  library(tibble)
})

# ---------------------------------------------------------------------
# Hallmark gene sets (frozen)
# ---------------------------------------------------------------------

if (!file.exists(file_hallmark)) {
  stop(
    "Hallmark gene sets not found at ", file_hallmark,
    "\nRun R/make_hallmark_sets.R once to create it.",
    call. = FALSE
  )
}

hallmark_sets <- read.delim(file_hallmark, stringsAsFactors = FALSE)
stopifnot(length(unique(hallmark_sets$gs_name)) == 50)

# ---------------------------------------------------------------------
# Build a ranked gene vector
#
# Duplicate gene symbols are collapsed by keeping the entry with the
# largest absolute value of the ranking metric.
# ---------------------------------------------------------------------

make_ranked_vector <- function(df, metric) {

  if (!metric %in% colnames(df)) {
    stop("Column '", metric, "' not found in input.", call. = FALSE)
  }

  r <- df %>%
    filter(!is.na(symbol), symbol != "", !is.na(.data[[metric]])) %>%
    group_by(symbol) %>%
    slice_max(order_by = abs(.data[[metric]]), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(.data[[metric]]))

  g <- r[[metric]]
  names(g) <- r$symbol
  sort(g, decreasing = TRUE)
}

# ---------------------------------------------------------------------
# Run GSEA once
# ---------------------------------------------------------------------

run_gsea <- function(df, metric, seed) {

  g <- make_ranked_vector(df, metric)

  set.seed(seed)

  x <- GSEA(
    geneList      = g,
    TERM2GENE     = hallmark_sets,
    pvalueCutoff  = 1,          # keep all pathways; filter downstream
    pAdjustMethod = "BH",
    minGSSize     = gsea_min_size,
    maxGSSize     = gsea_max_size,
    verbose       = FALSE,
    seed          = TRUE
  )

  res <- as.data.frame(x@result)
  res$seed        <- seed
  res$rank_metric <- metric
  res$n_ranked    <- length(g)
  res
}

# ---------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------

comparisons <- list(
  CRISPRa = file_crispra_de,
  DMC     = file_dmc_de
)

all_results <- list()

for (cmp in names(comparisons)) {

  path <- comparisons[[cmp]]

  if (!file.exists(path)) {
    stop("Input not found: ", path, call. = FALSE)
  }

  de <- read.csv(path, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)

  for (metric in gsea_rank_metrics) {

    message("GSEA: ", cmp, " ranked by ", metric)

    res <- bind_rows(lapply(gsea_seeds, function(s) run_gsea(de, metric, s)))

    tag <- paste0(cmp, "_rank", metric)
    all_results[[tag]] <- res

    write_csv_safe(
      res,
      file.path(results_dir, paste0("GSEA_Hallmark_", tag, "_allseeds.csv"))
    )

    primary <- res[res$seed == gsea_seeds[1], ]

    write_csv_safe(
      primary,
      file.path(results_dir, paste0("GSEA_Hallmark_", tag, ".csv"))
    )

    message(
      "  ", nrow(primary), " pathways; ",
      sum(primary$p.adjust < padj_cutoff), " with padj < ", padj_cutoff,
      "; ", unique(res$n_ranked), " ranked genes"
    )
  }
}

# ---------------------------------------------------------------------
# Sensitivity analysis: does the ranking metric change the conclusions?
#
# Compares NES and adjusted p-values between log2FoldChange-ranked and
# Wald-statistic-ranked runs at the primary seed, for each comparison.
# ---------------------------------------------------------------------

metric_comparison <- lapply(names(comparisons), function(cmp) {

  lfc  <- all_results[[paste0(cmp, "_ranklog2FoldChange")]]
  wald <- all_results[[paste0(cmp, "_rankstat")]]

  lfc  <- lfc[lfc$seed == gsea_seeds[1], c("ID", "NES", "p.adjust")]
  wald <- wald[wald$seed == gsea_seeds[1], c("ID", "NES", "p.adjust")]

  merged <- merge(lfc, wald, by = "ID", suffixes = c("_log2FC", "_Wald"))
  merged$comparison <- cmp

  merged$sig_log2FC <- merged$p.adjust_log2FC < padj_cutoff
  merged$sig_Wald   <- merged$p.adjust_Wald   < padj_cutoff
  merged$same_call  <- merged$sig_log2FC == merged$sig_Wald
  merged$same_sign  <- sign(merged$NES_log2FC) == sign(merged$NES_Wald)

  merged
}) %>% bind_rows()

write_csv_safe(
  metric_comparison,
  file.path(results_dir, "GSEA_rank_metric_comparison.csv")
)

for (cmp in unique(metric_comparison$comparison)) {
  m <- metric_comparison[metric_comparison$comparison == cmp, ]
  message(
    cmp, ": ", sum(m$same_call), "/", nrow(m),
    " pathways receive the same significance call, ",
    sum(m$same_sign), "/", nrow(m), " the same NES sign"
  )
}

save_session_info("02_gsea_hallmark")

message("02_gsea_hallmark.R complete.")
