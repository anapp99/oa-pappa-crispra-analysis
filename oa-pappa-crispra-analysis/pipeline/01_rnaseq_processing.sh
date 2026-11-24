#!/usr/bin/env bash
# =====================================================================
# 01_rnaseq_processing.sh
#
# Per-sample RNA-seq processing: quality assessment, trimming,
# alignment for QC, and transcript quantification.
#
# Usage:
#   ./01_rnaseq_processing.sh SAMPLE_NAME R1.fastq.gz R2.fastq.gz
#
# Or over all samples in a directory:
#   for r1 in fastq/*_R1.fastq.gz; do
#     s=$(basename "$r1" _R1.fastq.gz)
#     ./01_rnaseq_processing.sh "$s" "$r1" "${r1/_R1/_R2}"
#   done
#
# -------------------------------------------------------------------
# PROVENANCE NOTE
# This script documents the processing described in the Methods.
# Tool versions are as stated there. Lines marked [INFERRED] use tool
# defaults or conventional settings that the Methods do not specify;
# these were not necessarily the values used in the original run and
# should be confirmed before being described as such.
#
# Tools and versions from the Methods:
#   FastQC 0.12.1, MultiQC 1.32, fastp 0.23.4, STAR 2.7.10b,
#   SAMtools 1.18, Salmon 1.10.2, Picard, Preseq 2.0.1,
#   featureCounts 2.1.1
# -------------------------------------------------------------------
# =====================================================================

set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 SAMPLE_NAME R1.fastq.gz R2.fastq.gz" >&2
  exit 1
fi

SAMPLE="$1"
R1="$2"
R2="$3"

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

REF_DIR="${REF_DIR:-references}"
OUT_DIR="${OUT_DIR:-processed}"
THREADS="${THREADS:-8}"                      # [INFERRED] not specified

GENCODE_RELEASE=29
STAR_INDEX="${REF_DIR}/star_index_GRCh38_v${GENCODE_RELEASE}"
SALMON_INDEX="${REF_DIR}/salmon_index_gencode_v${GENCODE_RELEASE}_decoy"
GTF="${REF_DIR}/gencode.v${GENCODE_RELEASE}.annotation.gtf"

# NOTE: the Methods describe featureCounts using GENCODE v39 while
# Salmon used v29. If that is correct, point GTF_FEATURECOUNTS at the
# v39 annotation. Mixing annotation releases is worth confirming.
GTF_FEATURECOUNTS="${GTF_FEATURECOUNTS:-${GTF}}"

mkdir -p \
  "${OUT_DIR}/fastqc_raw" \
  "${OUT_DIR}/fastqc_trimmed" \
  "${OUT_DIR}/fastp" \
  "${OUT_DIR}/trimmed" \
  "${OUT_DIR}/star" \
  "${OUT_DIR}/salmon" \
  "${OUT_DIR}/qc"

echo "=== ${SAMPLE} ==="

# ---------------------------------------------------------------------
# 1. Quality assessment, raw reads
# ---------------------------------------------------------------------

fastqc \
  --threads "${THREADS}" \
  --outdir "${OUT_DIR}/fastqc_raw" \
  "${R1}" "${R2}"

# ---------------------------------------------------------------------
# 2. Adapter and quality trimming
#
# Automatic paired-end adapter detection; first 10 bases of each read
# trimmed to reduce sequence composition bias at read starts.
# Quality thresholds are fastp defaults. [INFERRED]
# ---------------------------------------------------------------------

T1="${OUT_DIR}/trimmed/${SAMPLE}_R1.trimmed.fastq.gz"
T2="${OUT_DIR}/trimmed/${SAMPLE}_R2.trimmed.fastq.gz"

fastp \
  --in1 "${R1}" \
  --in2 "${R2}" \
  --out1 "${T1}" \
  --out2 "${T2}" \
  --detect_adapter_for_pe \
  --trim_front1 10 \
  --trim_front2 10 \
  --thread "${THREADS}" \
  --json "${OUT_DIR}/fastp/${SAMPLE}.fastp.json" \
  --html "${OUT_DIR}/fastp/${SAMPLE}.fastp.html"

# ---------------------------------------------------------------------
# 3. Quality assessment, trimmed reads
# ---------------------------------------------------------------------

fastqc \
  --threads "${THREADS}" \
  --outdir "${OUT_DIR}/fastqc_trimmed" \
  "${T1}" "${T2}"

# ---------------------------------------------------------------------
# 4. Alignment (for QC metrics only)
#
# Quantification is performed independently by Salmon in step 5.
# These alignments support Picard and featureCounts metrics.
# ---------------------------------------------------------------------

STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX}" \
  --readFilesIn "${T1}" "${T2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUT_DIR}/star/${SAMPLE}." \
  --outSAMtype BAM Unsorted            # [INFERRED] sorted below with SAMtools

BAM_UNSORTED="${OUT_DIR}/star/${SAMPLE}.Aligned.out.bam"
BAM="${OUT_DIR}/star/${SAMPLE}.sorted.bam"

samtools sort -@ "${THREADS}" -o "${BAM}" "${BAM_UNSORTED}"
samtools index -@ "${THREADS}" "${BAM}"
rm -f "${BAM_UNSORTED}"

# ---------------------------------------------------------------------
# 5. Transcript quantification
#
# Run from the trimmed FASTQ files, independently of the STAR
# alignment, against the decoy-aware GENCODE v29 index.
#
# Library type is set to automatic detection. Bias correction flags
# (--gcBias, --seqBias) are NOT enabled here because the Methods do
# not mention them; enabling them changes quantification. [INFERRED]
# ---------------------------------------------------------------------

salmon quant \
  --index "${SALMON_INDEX}" \
  --libType A \
  --mates1 "${T1}" \
  --mates2 "${T2}" \
  --threads "${THREADS}" \
  --validateMappings \
  --output "${OUT_DIR}/salmon/${SAMPLE}"

# ---------------------------------------------------------------------
# 6. Post-alignment quality control
# ---------------------------------------------------------------------

# Duplicate rate
picard MarkDuplicates \
  I="${BAM}" \
  O="${OUT_DIR}/qc/${SAMPLE}.markdup.bam" \
  M="${OUT_DIR}/qc/${SAMPLE}.markdup_metrics.txt" \
  REMOVE_DUPLICATES=false              # [INFERRED] metrics only

# Library composition: coding, UTR, intronic, intergenic
# Requires a refFlat file derived from the annotation; build once with
# UCSC gtfToGenePred if not already present.
REFFLAT="${REF_DIR}/gencode.v${GENCODE_RELEASE}.refFlat.txt"

if [[ -f "${REFFLAT}" ]]; then
  picard CollectRnaSeqMetrics \
    I="${BAM}" \
    O="${OUT_DIR}/qc/${SAMPLE}.rnaseq_metrics.txt" \
    REF_FLAT="${REFFLAT}" \
    STRAND=NONE                        # [INFERRED] strandedness not stated
else
  echo "refFlat not found, skipping CollectRnaSeqMetrics: ${REFFLAT}" >&2
fi

# Library complexity
preseq lc_extrap \
  -bam \
  -output "${OUT_DIR}/qc/${SAMPLE}.preseq.txt" \
  "${BAM}" \
  || echo "preseq failed for ${SAMPLE}, continuing" >&2

# Gene assignment efficiency
featureCounts \
  -T "${THREADS}" \
  -p \
  -a "${GTF_FEATURECOUNTS}" \
  -o "${OUT_DIR}/qc/${SAMPLE}.featureCounts.txt" \
  "${BAM}"

echo "=== ${SAMPLE} complete ==="

# ---------------------------------------------------------------------
# After all samples, summarise with MultiQC:
#
#   multiqc "${OUT_DIR}" --outdir "${OUT_DIR}/multiqc"
#
# Then import the Salmon quantifications into R:
#
#   SALMON_DIR="${OUT_DIR}/salmon" Rscript R/01_deseq2_crispra.R
# ---------------------------------------------------------------------
