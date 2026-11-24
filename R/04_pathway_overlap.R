# =====================================================================
# 04_pathway_overlap.R
#
# Pathway-level comparison between CRISPRa vs SCR and DMC vs PLC.
# Produces the side-by-side Hallmark table underlying Figure 6D
# (Supplementary Table S9).
#
# Every Hallmark pathway is retained, with NES and adjusted p-values
# from both comparisons and flags for whether it is significant in
# both and whether the direction of enrichment agrees. Filtering is
# left to the reader rather than applied before export.
#
# Requires 02_gsea_hallmark.R to have been run first.
#
# Output: results/pathway_overlap_all_hallmark.csv   (all 50)
#         results/pathway_overlap_sig_both.csv       (significant in both)
#
# Run from the repository root:  Rscript R/04_pathway_overlap.R
# =====================================================================

source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ---------------------------------------------------------------------
# Read the primary-seed GSEA results for both comparisons
# ---------------------------------------------------------------------

gsea_path <- function(cmp) {
  file.path(
    results_dir,
    paste0("GSEA_Hallmark_", cmp, "_rank", gsea_primary_metric, ".csv")
  )
}

read_gsea <- function(cmp, suffix) {

  path <- gsea_path(cmp)

  if (!file.exists(path)) {
    stop(
      "GSEA results not found: ", path,
      "\nRun R/02_gsea_hallmark.R first.",
      call. = FALSE
    )
  }

  read_csv(path, show_col_types = FALSE) %>%
    select(ID, Description, NES, p.adjust, core_enrichment) %>%
    rename_with(
      ~ paste0(c("NES", "padj", "core"), "_", suffix),
      c("NES", "p.adjust", "core_enrichment")
    )
}

crispra <- read_gsea("CRISPRa", "CRISPRa")
dmc     <- read_gsea("DMC",     "DMC")

# ---------------------------------------------------------------------
# Join and annotate
# ---------------------------------------------------------------------

overlap <- full_join(crispra, dmc, by = c("ID", "Description")) %>%
  mutate(
    sig_CRISPRa = !is.na(padj_CRISPRa) & padj_CRISPRa < padj_cutoff,
    sig_DMC     = !is.na(padj_DMC)     & padj_DMC     < padj_cutoff,
    sig_both    = sig_CRISPRa & sig_DMC,
    direction = case_when(
      !sig_both                            ~ NA_character_,
      NES_CRISPRa > 0 & NES_DMC > 0        ~ "Up in both",
      NES_CRISPRa < 0 & NES_DMC < 0        ~ "Down in both",
      TRUE                                 ~ "Opposite direction"
    ),
    concordant = !is.na(direction) & direction != "Opposite direction",
    mean_NES   = (NES_CRISPRa + NES_DMC) / 2
  ) %>%
  arrange(desc(sig_both), desc(concordant), desc(mean_NES))

# ---------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------

# Complete table: all Hallmark pathways, nothing filtered out.
write_csv_safe(
  overlap,
  file.path(results_dir, "pathway_overlap_all_hallmark.csv")
)

sig_both <- filter(overlap, sig_both)

write_csv_safe(
  sig_both,
  file.path(results_dir, "pathway_overlap_sig_both.csv")
)

# ---------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------

message("\n", strrep("=", 60))
message("Hallmark pathways tested: ", nrow(overlap))
message("Significant in CRISPRa:   ", sum(overlap$sig_CRISPRa))
message("Significant in DMC vs PLC: ", sum(overlap$sig_DMC))
message("Significant in both:      ", nrow(sig_both))
message(strrep("-", 60))

if (nrow(sig_both) > 0) {
  print(
    sig_both %>%
      count(direction, name = "n_pathways")
  )

  discordant <- filter(sig_both, !concordant)

  if (nrow(discordant) > 0) {
    message("\nSignificant in both but opposite direction:")
    print(
      discordant %>%
        select(ID, NES_CRISPRa, padj_CRISPRa, NES_DMC, padj_DMC)
    )
  }
}

message(strrep("=", 60))

save_session_info("04_pathway_overlap")

message("04_pathway_overlap.R complete.")
