# =====================================================================
# 03_concordance_genes.R
#
# Gene-level overlap between CRISPRa vs SCR and DMC vs PLC cartilage.
#
# Two filters are applied:
#   (a) significance only        - padj < 0.05 in both comparisons.
#                                  This set underlies Figure 6C.
#   (b) significance + log2FC    - additionally |log2FC| > cutoff in
#                                  both, reported as a sensitivity check.
#
# Output: results/concordance_sig_only_all_shared.csv
#         results/concordance_sig_only_concordant.csv
#         results/concordance_sig_lfc_all_shared.csv
#         results/concordance_sig_lfc_concordant.csv
#         results/concordance_summary.csv
#
# Run from the repository root:  Rscript R/03_concordance_genes.R
# =====================================================================

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

required_cols <- c("symbol", "log2FoldChange", "padj")

# ---------------------------------------------------------------------
# Read and standardise one differential expression table
# ---------------------------------------------------------------------

load_de <- function(path, suffix) {

  if (!file.exists(path)) {
    stop("Input not found: ", path, call. = FALSE)
  }

  df <- read_csv(path, show_col_types = FALSE)

  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(
      "Missing columns in ", basename(path), ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  df %>%
    select(all_of(required_cols)) %>%
    filter(!is.na(symbol), symbol != "") %>%
    # Collapse duplicate symbols, keeping the most significant entry
    arrange(padj) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    rename_with(
      ~ paste0(.x, "_", suffix),
      c("log2FoldChange", "padj")
    )
}

# ---------------------------------------------------------------------
# Overlap and concordance
# ---------------------------------------------------------------------

run_concordance <- function(exp_df, oa_df, label, lfc_cutoff = NULL) {

  message("\n", strrep("=", 60))
  message(label)
  message(
    if (is.null(lfc_cutoff)) {
      "Filter: padj only"
    } else {
      paste0("Filter: padj and |log2FC| > ", lfc_cutoff)
    }
  )
  message(strrep("=", 60))

  exp_sig <- exp_df %>% filter(!is.na(padj_EXP), padj_EXP < padj_cutoff)
  oa_sig  <- oa_df  %>% filter(!is.na(padj_OA),  padj_OA  < padj_cutoff)

  if (!is.null(lfc_cutoff)) {
    exp_sig <- exp_sig %>% filter(abs(log2FoldChange_EXP) > lfc_cutoff)
    oa_sig  <- oa_sig  %>% filter(abs(log2FoldChange_OA)  > lfc_cutoff)
  }

  message("CRISPRa significant genes: ", nrow(exp_sig))
  message("DMC vs PLC significant genes: ", nrow(oa_sig))

  shared <- inner_join(exp_sig, oa_sig, by = "symbol") %>%
    mutate(
      direction_EXP = case_when(
        log2FoldChange_EXP > 0 ~ "Up",
        log2FoldChange_EXP < 0 ~ "Down",
        TRUE                   ~ "No change"
      ),
      direction_OA = case_when(
        log2FoldChange_OA > 0 ~ "Up",
        log2FoldChange_OA < 0 ~ "Down",
        TRUE                  ~ "No change"
      ),
      concordant = direction_EXP == direction_OA,
      mean_log2FC = (log2FoldChange_EXP + log2FoldChange_OA) / 2
    ) %>%
    arrange(desc(concordant), desc(abs(mean_log2FC)))

  n_shared     <- nrow(shared)
  n_up_both    <- sum(shared$direction_EXP == "Up"   & shared$direction_OA == "Up")
  n_down_both  <- sum(shared$direction_EXP == "Down" & shared$direction_OA == "Down")
  n_concordant <- sum(shared$concordant)
  pct          <- if (n_shared > 0) round(100 * n_concordant / n_shared, 1) else NA

  message("Significant in both: ", n_shared)
  message("Up in both: ",   n_up_both)
  message("Down in both: ", n_down_both)
  message("Concordant: ",   n_concordant, " (", pct, "%)")

  list(
    shared = shared,
    summary = data.frame(
      analysis       = label,
      lfc_cutoff     = if (is.null(lfc_cutoff)) NA_real_ else lfc_cutoff,
      n_exp_sig      = nrow(exp_sig),
      n_oa_sig       = nrow(oa_sig),
      n_shared       = n_shared,
      n_up_both      = n_up_both,
      n_down_both    = n_down_both,
      n_concordant   = n_concordant,
      pct_concordant = pct
    )
  )
}

# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------

exp_df <- load_de(file_crispra_de, "EXP")
oa_df  <- load_de(file_dmc_de,     "OA")

runs <- list(
  sig_only = run_concordance(
    exp_df, oa_df,
    label = "CRISPRa vs DMC/PLC - significance only",
    lfc_cutoff = NULL
  ),
  sig_lfc = run_concordance(
    exp_df, oa_df,
    label = "CRISPRa vs DMC/PLC - significance and fold change",
    lfc_cutoff = 1
  )
)

for (tag in names(runs)) {

  shared <- runs[[tag]]$shared

  write_csv_safe(
    shared,
    file.path(results_dir, paste0("concordance_", tag, "_all_shared.csv"))
  )

  write_csv_safe(
    filter(shared, concordant),
    file.path(results_dir, paste0("concordance_", tag, "_concordant.csv"))
  )
}

summary_all <- bind_rows(lapply(runs, function(x) x$summary))

write_csv_safe(
  summary_all,
  file.path(results_dir, "concordance_summary.csv")
)

message("\nSummary:")
print(summary_all)

save_session_info("03_concordance_genes")

message("03_concordance_genes.R complete.")
