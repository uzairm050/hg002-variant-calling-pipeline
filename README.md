# HG002 Variant Calling Pipeline

## Project Summary

This project implements a complete DNA variant calling pipeline on the HPC cluster using containerized bioinformatics tools. We process PacBio HiFi sequencing reads from HG002 (a reference human sample maintained by GIAB), identify genetic variants using two independent tools, and evaluate their accuracy against a validated truth set.

A genetic variant is a position in the genome where an individual's DNA differs from the reference human genome. Identifying these variants accurately is fundamental to genomic research and clinical diagnostics.

## Pipeline Overview

The pipeline consists of five stages:

1. **Read Alignment** — Raw sequencing reads are aligned to the GRCh38 human reference genome using Minimap2, which is optimized for long PacBio HiFi reads.

2. **BAM Processing** — The alignment file is sorted and indexed using Samtools to prepare it for variant calling.

3. **Variant Calling with Clair3** — A fast, neural network-based variant caller that uses a two-stage approach: pileup calling followed by full-alignment calling for low-confidence sites.

4. **Variant Calling with DeepVariant** — Google's deep learning variant caller that treats each candidate variant as an image classification problem using a convolutional neural network.

5. **Benchmarking with hap.py** — Both variant call sets are evaluated against the GIAB HG002 v4.2.1 truth set to measure precision, recall, and F1 score.

## How to Run

The entire pipeline is packaged as a SLURM job script. To run it on the HPC cluster:
sbatch pipeline.sh
SLURM will allocate the required compute resources (8 CPUs, 32GB RAM) and execute all steps automatically.

## Results

Note: This analysis used 25% of the full HG002 dataset. Low recall is therefore expected, as many genomic regions lacked sufficient read coverage for variant detection. At full coverage, both tools typically achieve over 99% recall.

### SNP Performance

| Metric | Clair3 | DeepVariant |
|--------|--------|-------------|
| Variants Called | 255,292 | 181,989 |
| True Positives | 167,520 | 101,556 |
| Precision | 82.5% | 76.2% |
| Recall | 4.98% | 3.02% |
| F1 Score | 0.094 | 0.058 |

### INDEL Performance

| Metric | Clair3 | DeepVariant |
|--------|--------|-------------|
| Variants Called | 39,154 | 33,105 |
| True Positives | 14,324 | 10,792 |
| Precision | 62.6% | 64.9% |
| Recall | 2.73% | 2.05% |
| F1 Score | 0.052 | 0.040 |

### Key Finding

Clair3 outperformed DeepVariant at this coverage level, achieving higher recall and precision for SNPs. DeepVariant showed marginally better INDEL precision. At full coverage, results are expected to be more comparable.

## Tools and Containers

| Tool | Version | Container | Purpose |
|------|---------|-----------|---------|
| Minimap2 | latest | docker://staphb/minimap2 | Read alignment |
| Samtools | latest | docker://staphb/samtools | BAM processing |
| Clair3 | latest | docker://hkubal/clair3 | Variant calling |
| DeepVariant | 1.6.1 | docker://google/deepvariant | Variant calling |
| hap.py | latest | docker://pkrusche/hap.py | Benchmarking |

## Repository Contents

| File | Description |
|------|-------------|
| pipeline.sh | SLURM job script — runs the full pipeline end to end |
| deepvariant.nf | Nextflow workflow script for DeepVariant |
