#!/usr/bin/env bash
set -euo pipefail

# HG002 Variant Calling Pipeline
# Workflow: minimap2 (map-hifi) -> samtools sort/index -> Clair3 + DeepVariant -> VCF -> Benchmarking
#
# Usage:
#   bash runpipeline.sh                        # uses current directory as BASEDIR
#   bash runpipeline.sh /path/to/workdir       # uses specified directory as BASEDIR
#
# Expected directory structure:
#   $BASEDIR/
#   ├── data/
#   │   ├── GRCh38.fa                                        # GRCh38 reference genome
#   │   ├── HG002_quarter.fastq.gz                           # HG002 HiFi quarter subset
#   │   └── giab/
#   │       ├── HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz   # GIAB truth VCF
#   │       └── HG002_GRCh38_1_22_v4.2.1_benchmark.bed      # GIAB confident regions BED
#   └── models/
#       └── hifi_revio/                                      # Clair3 hifi_revio model directory
#
# GIAB truth files can be downloaded from:
#   https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/

BASEDIR=${1:-$(pwd)}

# Validate expected inputs exist before launching
if [ ! -f "${BASEDIR}/data/GRCh38.fa" ]; then
  echo "ERROR: Reference not found at ${BASEDIR}/data/GRCh38.fa" >&2
  exit 1
fi
if [ ! -f "${BASEDIR}/data/HG002_quarter.fastq.gz" ]; then
  echo "ERROR: Reads not found at ${BASEDIR}/data/HG002_quarter.fastq.gz" >&2
  exit 1
fi
if [ ! -d "${BASEDIR}/models/hifi_revio" ]; then
  echo "ERROR: Clair3 model not found at ${BASEDIR}/models/hifi_revio" >&2
  exit 1
fi

echo "======================================"
echo " HG002 Variant Calling Pipeline"
echo " BASEDIR : ${BASEDIR}"
echo " Step 1/4: Running Clair3 pipeline"
echo "======================================"

nextflow run main.nf -profile slurm \
  --reference  "${BASEDIR}/data/GRCh38.fa" \
  --reads      "${BASEDIR}/data/HG002_quarter.fastq.gz" \
  --model_path "${BASEDIR}/models/hifi_revio" \
  --sample     HG002 \
  --outdir     "${BASEDIR}/results/clair3" \
  --threads    8

echo ""
echo "======================================"
echo " Step 2/4: Running DeepVariant pipeline"
echo "======================================"

nextflow run deepvariant.nf -profile slurm \
  --reference  "${BASEDIR}/data/GRCh38.fa" \
  --reads      "${BASEDIR}/data/HG002_quarter.fastq.gz" \
  --bam        "${BASEDIR}/results/clair3/HG002.sorted.bam" \
  --sample     HG002 \
  --outdir     "${BASEDIR}/results/deepvariant" \
  --threads    8

echo ""
echo "======================================"
echo " Step 3/4: BCFtools — filter to chr1-22"
echo "======================================"

CLAIR3_VCF="${BASEDIR}/results/clair3/HG002.vcf.gz"
DV_VCF="${BASEDIR}/results/deepvariant/HG002.vcf.gz"
CLAIR3_FILTERED="${BASEDIR}/results/clair3/HG002.chr1-22.vcf.gz"
DV_FILTERED="${BASEDIR}/results/deepvariant/HG002.chr1-22.vcf.gz"

# Build chr1-22 region list
REGIONS=$(seq 1 22 | sed 's/^/chr/' | tr '\n' ',' | sed 's/,$//')

singularity exec docker://staphb/bcftools:latest \
  bcftools view \
    --regions "${REGIONS}" \
    --output-type z \
    --output "${CLAIR3_FILTERED}" \
    "${CLAIR3_VCF}"

singularity exec docker://staphb/bcftools:latest \
  bcftools index --tbi "${CLAIR3_FILTERED}"

echo "Clair3 filtered VCF: ${CLAIR3_FILTERED}"

singularity exec docker://staphb/bcftools:latest \
  bcftools view \
    --regions "${REGIONS}" \
    --output-type z \
    --output "${DV_FILTERED}" \
    "${DV_VCF}"

singularity exec docker://staphb/bcftools:latest \
  bcftools index --tbi "${DV_FILTERED}"

echo "DeepVariant filtered VCF: ${DV_FILTERED}"

echo ""
echo "======================================"
echo " Step 4/4: hap.py benchmarking vs GIAB NISTv4.2.1"
echo "======================================"

# Truth set — expected at $BASEDIR/data/giab/
TRUTH_VCF="${BASEDIR}/data/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${BASEDIR}/data/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.bed"
REF="${BASEDIR}/data/GRCh38.fa"

# Validate truth set exists
if [ ! -f "${TRUTH_VCF}" ]; then
  echo "ERROR: GIAB truth VCF not found at ${TRUTH_VCF}" >&2
  echo "Download from: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/" >&2
  exit 1
fi
if [ ! -f "${TRUTH_BED}" ]; then
  echo "ERROR: GIAB truth BED not found at ${TRUTH_BED}" >&2
  exit 1
fi

mkdir -p "${BASEDIR}/results/benchmark/clair3"
mkdir -p "${BASEDIR}/results/benchmark/deepvariant"

echo "Benchmarking Clair3..."
singularity exec docker://pkrusche/hap.py:latest \
  /opt/hap.py/bin/hap.py \
    "${TRUTH_VCF}" \
    "${CLAIR3_FILTERED}" \
    --reference   "${REF}" \
    --confident-regions "${TRUTH_BED}" \
    --engine      vcfeval \
    --threads     8 \
    --output-file "${BASEDIR}/results/benchmark/clair3/HG002_clair3"

echo "Benchmarking DeepVariant..."
singularity exec docker://pkrusche/hap.py:latest \
  /opt/hap.py/bin/hap.py \
    "${TRUTH_VCF}" \
    "${DV_FILTERED}" \
    --reference   "${REF}" \
    --confident-regions "${TRUTH_BED}" \
    --engine      vcfeval \
    --threads     8 \
    --output-file "${BASEDIR}/results/benchmark/deepvariant/HG002_deepvariant"

echo ""
echo "======================================"
echo " ALL STEPS COMPLETE"
echo "======================================"
echo " Clair3 VCF (filtered)    : ${CLAIR3_FILTERED}"
echo " DeepVariant VCF (filtered): ${DV_FILTERED}"
echo " Clair3 benchmark         : ${BASEDIR}/results/benchmark/clair3/HG002_clair3.summary.csv"
echo " DeepVariant benchmark    : ${BASEDIR}/results/benchmark/deepvariant/HG002_deepvariant.summary.csv"
echo ""
echo " Results directory layout:"
echo "   results/"
echo "   ├── clair3/             Clair3 BAM + VCF"
echo "   ├── deepvariant/        DeepVariant VCF"
echo "   └── benchmark/"
echo "       ├── clair3/         hap.py precision/recall vs GIAB"
echo "       └── deepvariant/    hap.py precision/recall vs GIAB"
echo "======================================"
