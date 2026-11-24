#!/usr/bin/env bash
# =====================================================================
# 00_build_references.sh
#
# Builds the STAR genome index and the decoy-aware Salmon index used by
# 01_rnaseq_processing.sh.
#
# Reference: GRCh38 primary assembly, GENCODE v29 transcriptome.
#
# -------------------------------------------------------------------
# PROVENANCE NOTE
# This script documents the reference build implied by the Methods.
# Lines marked [INFERRED] use tool defaults or conventional settings
# that the Methods do not specify. Confirm against the original run
# before describing these as the parameters used.
# -------------------------------------------------------------------
# =====================================================================

set -euo pipefail

# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

REF_DIR="${REF_DIR:-references}"
THREADS="${THREADS:-8}"                      # [INFERRED] not specified
READ_LENGTH="${READ_LENGTH:-150}"            # 150 bp paired-end

GENCODE_RELEASE=29
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"

GENOME_FA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
TX_FA="${REF_DIR}/gencode.v${GENCODE_RELEASE}.transcripts.fa"
GTF="${REF_DIR}/gencode.v${GENCODE_RELEASE}.annotation.gtf"

STAR_INDEX="${REF_DIR}/star_index_GRCh38_v${GENCODE_RELEASE}"
SALMON_INDEX="${REF_DIR}/salmon_index_gencode_v${GENCODE_RELEASE}_decoy"

mkdir -p "${REF_DIR}"

# ---------------------------------------------------------------------
# 1. Download references
# ---------------------------------------------------------------------

download() {
  local url="$1" out="$2"
  if [[ -f "${out}" ]]; then
    echo "Already present, skipping: ${out}"
    return
  fi
  echo "Downloading: ${url}"
  curl -fsSL "${url}" -o "${out}.gz"
  gunzip "${out}.gz"
}

download "${GENCODE_BASE}/GRCh38.primary_assembly.genome.fa.gz"        "${GENOME_FA}"
download "${GENCODE_BASE}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz" "${TX_FA}"
download "${GENCODE_BASE}/gencode.v${GENCODE_RELEASE}.annotation.gtf.gz" "${GTF}"

# ---------------------------------------------------------------------
# 2. STAR index
#
# Used for alignment-based QC metrics only (Picard, featureCounts).
# Quantification is done independently by Salmon.
# ---------------------------------------------------------------------

if [[ ! -d "${STAR_INDEX}" ]]; then
  echo "Building STAR index..."
  mkdir -p "${STAR_INDEX}"
  STAR \
    --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --genomeFastaFiles "${GENOME_FA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang $((READ_LENGTH - 1))   # [INFERRED] conventional value
else
  echo "STAR index already present, skipping."
fi

# ---------------------------------------------------------------------
# 3. Salmon decoy-aware index
#
# Decoys are the full genome, following the Salmon documentation's
# recommended approach for a decoy-aware index.
# ---------------------------------------------------------------------

if [[ ! -d "${SALMON_INDEX}" ]]; then
  echo "Building Salmon decoy-aware index..."

  # Decoy sequence names = genome contig names
  grep "^>" "${GENOME_FA}" | cut -d " " -f 1 | sed 's/>//g' \
    > "${REF_DIR}/decoys.txt"

  # Transcriptome followed by genome
  cat "${TX_FA}" "${GENOME_FA}" > "${REF_DIR}/gentrome.fa"

  salmon index \
    --transcripts "${REF_DIR}/gentrome.fa" \
    --decoys "${REF_DIR}/decoys.txt" \
    --index "${SALMON_INDEX}" \
    --gencode \
    --threads "${THREADS}" \
    --kmerLen 31                          # [INFERRED] Salmon default
else
  echo "Salmon index already present, skipping."
fi

echo "Reference build complete."
echo "  STAR index:   ${STAR_INDEX}"
echo "  Salmon index: ${SALMON_INDEX}"
