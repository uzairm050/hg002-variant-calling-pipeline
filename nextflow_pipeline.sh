#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --output=logs/nextflow_%j.log
#SBATCH --error=logs/nextflow_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=16G

export PATH=/home/apps/nextflow_framework/framework/25.10.4:$PATH
export NXF_HOME=$(pwd)/.nextflow_home

echo "Starting pipeline at $(date)"
echo "Working directory: $(pwd)"

nextflow run scripts/variant_calling.nf \
  -with-singularity $(pwd)/deepvariant_1.6.1.sif

echo "Pipeline finished at $(date)"
