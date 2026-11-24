# =====================================================================
# make_hallmark_sets.R
#
# Writes a frozen copy of the MSigDB Hallmark gene sets to
# data/hallmark_sets.tsv, so that 02_gsea_hallmark.R does not depend
# on the installed msigdbr version. Run once; the resulting file is
# committed to the repository.
#
# Run from the repository root:  Rscript R/make_hallmark_sets.R
# =====================================================================

source("config.R")

suppressPackageStartupMessages({
  library(msigdbr)
  library(dplyr)
})

if (file.exists(file_hallmark)) {
  stop(
    "Hallmark file already exists: ", file_hallmark,
    "\nDelete it first if you intend to regenerate it.",
    call. = FALSE
  )
}

# msigdbr changed its argument name from `category` to `collection`
# in version 10. Support both.
hallmark <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "H")
)

hallmark_sets <- hallmark %>%
  select(gs_name, gene_symbol) %>%
  distinct() %>%
  arrange(gs_name, gene_symbol)

stopifnot(length(unique(hallmark_sets$gs_name)) == 50)

dir.create(dirname(file_hallmark), showWarnings = FALSE, recursive = TRUE)
write.table(
  hallmark_sets,
  file_hallmark,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message(
  "Wrote ", file_hallmark, ": ",
  length(unique(hallmark_sets$gs_name)), " gene sets, ",
  nrow(hallmark_sets), " gene-set memberships"
)

message("msigdbr version: ", as.character(packageVersion("msigdbr")))

save_session_info("make_hallmark_sets")
