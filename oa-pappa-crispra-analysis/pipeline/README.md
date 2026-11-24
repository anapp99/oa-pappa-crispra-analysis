# RNA-seq processing

Shell pipeline taking raw paired-end FASTQ files to the Salmon
quantifications used by `R/01_deseq2_crispra.R`.

## Status of this code

These scripts **document the processing described in the manuscript
Methods**. They were written from that description rather than exported
from the original run, so they reproduce the documented steps, tools and
versions, but not necessarily every parameter.

Settings the Methods do not specify are marked `[INFERRED]` in the
scripts and use tool defaults or conventional values. The ones most
likely to affect results:

| Setting | Here | Why it matters |
|---|---|---|
| Salmon `--gcBias` / `--seqBias` | Not enabled | Bias correction changes abundance estimates |
| Salmon `--libType` | `A` (auto-detect) | A fixed library type would differ |
| fastp quality thresholds | Defaults | Affects how much is trimmed |
| STAR `--sjdbOverhang` | `readLength - 1` | Affects splice junction detection |
| `CollectRnaSeqMetrics STRAND` | `NONE` | Strandedness is not stated in the Methods |

If you can recover the original invocations, replace the inferred values
and delete this section.

## Annotation release

Salmon quantifies against **GENCODE v29**, as stated in the Methods.
The Methods also describe featureCounts using **GENCODE v39** for QC.
`GTF_FEATURECOUNTS` defaults to v29 here; set it explicitly if v39 was
genuinely used, and note that mixing releases is worth confirming since
gene models differ between them.

## Running

```bash
# Once: download references and build indices
./pipeline/00_build_references.sh

# Per sample
./pipeline/01_rnaseq_processing.sh SAMPLE_NAME R1.fastq.gz R2.fastq.gz

# Or all samples
for r1 in fastq/*_R1.fastq.gz; do
  s=$(basename "$r1" _R1.fastq.gz)
  ./pipeline/01_rnaseq_processing.sh "$s" "$r1" "${r1/_R1/_R2}"
done

# Summarise QC
multiqc processed --outdir processed/multiqc

# Into R
SALMON_DIR=processed/salmon Rscript R/01_deseq2_crispra.R
```

Override `REF_DIR`, `OUT_DIR` and `THREADS` as environment variables.

## On a cluster

Each sample is independent, so the per-sample script suits an array job.
Adjust `THREADS` to match the cores requested.

## Software versions

As stated in the Methods: FastQC 0.12.1, MultiQC 1.32, fastp 0.23.4,
STAR 2.7.10b, SAMtools 1.18, Salmon 1.10.2, Preseq 2.0.1,
featureCounts 2.1.1, Picard. Downstream R analysis uses R 4.3.3 with
tximeta 1.20.3 and DESeq2 1.42.1.
