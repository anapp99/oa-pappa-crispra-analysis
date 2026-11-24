# Functional interrogation of osteoarthritis-associated regulatory variants

Analysis code for the RNA-sequencing and pathway analyses reported in
Piedade *et al.*, *Functional interrogation of osteoarthritis-associated
regulatory variants and effector genes*.

This repository covers the CRISPRa RNA-seq analysis and its comparison with
publicly available osteoarthritis cartilage transcriptomic data. Upstream
sequence processing pipelines are in separate repositories, listed below.

## What this produces

| Script | Produces |
|---|---|
| `pipeline/00_build_references.sh` | STAR and decoy-aware Salmon indices (GRCh38, GENCODE v29) |
| `pipeline/01_rnaseq_processing.sh` | FASTQ to Salmon quantifications, with QC |
| `R/01_deseq2_crispra.R` | Differential expression, CRISPRa vs SCR (Supplementary Table S6; Figure 6A; Supplementary Figure 4) |
| `R/02_gsea_hallmark.R` | Hallmark GSEA for both comparisons (Supplementary Tables S7, S8; Figure 6B) |
| `R/03_concordance_genes.R` | Gene-level overlap between comparisons (Figure 6C) |
| `R/04_pathway_overlap.R` | Pathway-level overlap (Supplementary Table S9; Figure 6D) |
| `R/make_hallmark_sets.R` | One-off: freezes the Hallmark gene sets to `data/` |

Sequence processing is documented in `pipeline/`. See `pipeline/README.md`
for the status of those scripts, which were written from the manuscript
Methods rather than exported from the original run.

## Running

Scripts are run from the repository root, in order:

```bash
Rscript R/make_hallmark_sets.R    # once only
Rscript R/01_deseq2_crispra.R
Rscript R/02_gsea_hallmark.R
Rscript R/03_concordance_genes.R
Rscript R/04_pathway_overlap.R
```

All paths and thresholds live in `config.R`; no other file contains a path.
Scripts 02 to 04 start from the differential expression tables in `data/`,
so they can be run without the raw sequencing data.

Script 01 additionally needs Salmon quantifications. Point it at them with:

```bash
export SALMON_DIR=/path/to/salmon/quants
```

Output files are written to `results/`. Scripts refuse to overwrite an
existing output, so delete `results/` to re-run cleanly.

## Inputs expected in `data/`

| File | Description |
|---|---|
| `PAPPA_RNAseq_Meta.csv` | Sample sheet with `Sample_Name`, `Case_Control` (C = SCR, E = CRISPRa), `Cell_Type` |
| `DE_DMC_vs_PLC_paired.csv` | Differential expression for damaged vs preserved cartilage, from published data |
| `hallmark_sets.tsv` | Written by `make_hallmark_sets.R` |

`DE_PAPPA_CRISPRa_vs_SCR.csv` is written by script 01 and read by 02 and 03.

## Analysis notes

**Ranking metric.** Pre-ranked GSEA is run twice: ranked by
`log2FoldChange`, which produced the results presented in the paper, and
ranked by the DESeq2 Wald statistic (`stat`) as a sensitivity analysis.
`results/GSEA_rank_metric_comparison.csv` reports how far the two agree.

**Seeds.** `clusterProfiler::GSEA` uses permutation testing, so results vary
slightly between runs. Each analysis is run under three seeds (123, 1, 42).
Seed 123 is primary; the other two confirm stability.

**Fold changes.** DESeq2 log2 fold changes are maximum-likelihood estimates
without shrinkage, matching the values quoted in the manuscript. These
differ slightly from the difference between group medians of normalised
counts, as plotted in Supplementary Figure 4.

**Gene sets.** Hallmark gene sets are read from a committed local file
rather than queried from MSigDB at run time, so results do not change when
MSigDB is updated.

**Complete results.** Supplementary tables are exported unfiltered. Pathway
and gene tables include non-significant and discordant entries, with flags,
rather than only the subset shown in the figures.

## Related repositories

| Pipeline | Repository |
|---|---|
| ATAC-seq processing | https://github.com/jakebutler98/UoM_ATAC_processor |
| CUT&Tag processing | https://github.com/CarlosRangel23/CutTagFlow |
| MPRA analysis | https://github.com/srstuppy/MPRA |

## Software

Analyses were run in R. Exact package versions for each script are recorded
in `results/sessionInfo_*.txt` when the scripts are run.

Required packages: DESeq2, tximeta, clusterProfiler, msigdbr, org.Hs.eg.db,
AnnotationDbi, dplyr, readr, tibble.

## Citation

If you use this code, please cite the paper. The archived release has a DOI
via Zenodo (see the repository sidebar).

## Licence

MIT. See `LICENSE`.
