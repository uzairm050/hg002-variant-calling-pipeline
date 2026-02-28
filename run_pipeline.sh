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
#   │   ├── GRCh38.fa                          # GRCh38 reference genome
#   │   └── HG002_quarter.fastq.gz             # HG002 HiFi quarter subset
#   └── models/
#       └── hifi_revio/                        # Clair3 hifi_revio model directory

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
echo " Step 1/2: Running Clair3 pipeline"
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
echo " Step 2/2: Running DeepVariant pipeline"
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
echo " Pipeline complete!"
echo " Clair3 VCF    : ${BASEDIR}/results/clair3/HG002.vcf.gz"
echo " DeepVariant VCF: ${BASEDIR}/results/deepvariant/HG002.vcf.gz"
echo " Next: run BCFtools filter + hap.py benchmarking"
echo "======================================"
