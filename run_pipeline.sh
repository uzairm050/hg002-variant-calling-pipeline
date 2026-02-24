#!/usr/bin/env bash
set -euo pipefail

# HG002 Variant Calling Pipeline
# Workflow: minimap2 (map-hifi) -> samtools sort/index -> Clair3 -> VCF
#
# Paths configured for this HPC cluster:
#   Reference  : GRCh38 human reference genome
#   Reads      : HG002 PacBio HiFi quarter subset
#   Model      : Clair3 hifi_revio model

BASEDIR=/hdd4/sines/advancedcomputationalbiology/uzair.sines/acbassignment

nextflow run main.nf -profile singularity \
  --reference ${BASEDIR}/data/GRCh38.fa \
  --reads ${BASEDIR}/data/HG002_quarter.fastq.gz \
  --model_path /opt/models/hifi_revio \
  --sample HG002 \
  --outdir ${BASEDIR}/results \
  --threads 8
