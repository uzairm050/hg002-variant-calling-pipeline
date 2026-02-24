#!/usr/bin/env bash
set -euo pipefail

# HG002 Variant Calling Pipeline
# Workflow: minimap2 (map-hifi) -> samtools sort/index -> Clair3 -> VCF
#
# Edit the paths below before running:
#   --reference  : path to your GRCh38 reference FASTA
#   --reads      : path to your HiFi FASTQ/FASTQ.GZ
#   --model_path : path to Clair3 HiFi model directory

nextflow run main.nf -profile singularity \
  --reference data/GRCh38.fa \
  --reads data/HG002_quarter.fastq.gz \
  --model_path /opt/models/hifi_revio \
  --sample HG002 \
  --outdir results \
  --threads 8
