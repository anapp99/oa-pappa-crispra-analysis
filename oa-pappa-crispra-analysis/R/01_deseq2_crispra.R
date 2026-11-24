# =====================================================================
# 01_deseq2_crispra.R
#
# Differential expression: CRISPRa of the rs1895062 enhancer versus
# non-targeting scrambled control (SCR), in SW1353 cells.
#
# Input : Salmon quantifications (one directory per sample) + metadata
# Output: data/DE_PAPPA_CRISPRa_vs_SCR.csv   (full results, all genes)
#         results/PAPPA_TAD_genes.csv        (PAPPA, ASTN2, TRIM32)
#         results/PAPPA_normalised_counts.csv
#         results/PAPPA_vst.rds, PAPPA_dds.rds
#
# Run from the repository root:  Rscript R/01_deseq2_crispra.R
# =====================================================================

source("config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximeta)
  library(readr)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ---------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------

meta <- read_csv(meta_file, show_col_types = FALSE)

# Case_Control: C = SCR control, E = CRISPRa
meta$Case_Control <- factor(meta$Case_Control, levels = c("C", "E"))

message("Sample groups:")
print(table(meta$Case_Control))

# ---------------------------------------------------------------------
# tximeta sample sheet
# ---------------------------------------------------------------------

coldata <- data.frame(
  names        = meta$Sample_Name,
  files        = file.path(salmon_dir, meta$Sample_Name, "quant.sf"),
  Case_Control = meta$Case_Control,
  Cell_Type    = meta$Cell_Type,
  stringsAsFactors = FALSE
)

missing <- coldata$files[!file.exists(coldata$files)]
if (length(missing) > 0) {
  stop(
    "Salmon quant.sf files not found:\n  ",
    paste(missing, collapse = "\n  "),
    "\nSet SALMON_DIR or edit salmon_dir in config.R.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# Import and summarise to gene level
# ---------------------------------------------------------------------

se  <- tximeta(coldata)
gse <- summarizeToGene(se)

# ---------------------------------------------------------------------
# DESeq2
# ---------------------------------------------------------------------

dds <- DESeqDataSet(gse, design = ~ Case_Control)

# Pre-filter lowly expressed genes
dds <- dds[rowMeans(counts(dds)) >= min_mean_count, ]
message("Genes retained after filtering: ", nrow(dds))

dds <- DESeq(dds)

# MLE log2 fold changes (no shrinkage), matching the values reported
# in the manuscript and figure legends.
res <- results(dds, contrast = c("Case_Control", "E", "C"))

# ---------------------------------------------------------------------
# Annotate
# ---------------------------------------------------------------------

res_df <- as.data.frame(res)
res_df$ensembl <- sub("\\..*", "", rownames(res_df))

res_df$symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = res_df$ensembl,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Fall back to the Ensembl ID where no symbol is available
no_symbol <- is.na(res_df$symbol)
res_df$symbol[no_symbol] <- res_df$ensembl[no_symbol]

res_df <- res_df %>%
  relocate(symbol, ensembl, .before = log2FoldChange) %>%
  arrange(padj)

# ---------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------

# Full results for every gene tested. This is Supplementary Table S6
# and the input to scripts 02 and 03.
write_csv_safe(res_df, file_crispra_de)

n_sig <- sum(
  res_df$padj < padj_cutoff &
    abs(res_df$log2FoldChange) > lfc_cutoff_volcano,
  na.rm = TRUE
)

message(
  "Significant genes (padj < ", padj_cutoff,
  " and |log2FC| > ", lfc_cutoff_volcano, "): ", n_sig
)

# Genes within the TAD containing rs1895062
tad_genes <- res_df %>%
  filter(symbol %in% c("PAPPA", "ASTN2", "TRIM32")) %>%
  select(symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

print(tad_genes)
write_csv_safe(tad_genes, file.path(results_dir, "PAPPA_TAD_genes.csv"))

# Normalised counts (Supplementary Figure 4)
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame()

norm_counts$ensembl <- sub("\\..*", "", rownames(norm_counts))
norm_counts$symbol <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys      = norm_counts$ensembl,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)
norm_counts <- relocate(norm_counts, symbol, ensembl)

write_csv_safe(
  norm_counts,
  file.path(results_dir, "PAPPA_normalised_counts.csv")
)

# Objects for downstream plotting
saveRDS(vst(dds, blind = FALSE), file.path(results_dir, "PAPPA_vst.rds"))
saveRDS(dds, file.path(results_dir, "PAPPA_dds.rds"))

save_session_info("01_deseq2_crispra")

message("01_deseq2_crispra.R complete.")
