# HG002 Variant Calling Pipeline

## Project Summary

This project implements a complete short variant calling pipeline on an HPC cluster using SLURM, Nextflow, and Singularity containers. We process PacBio HiFi sequencing reads from HG002 (a reference human sample maintained by GIAB), identify genetic variants across chromosomes 1-22 using two independent tools, and evaluate their accuracy against a validated truth set.

A genetic variant is a position in the genome where an individual's DNA differs from the reference human genome. Identifying these variants accurately is fundamental to genomic research and clinical diagnostics.

## Pipeline Architecture

The pipeline uses three core technologies working together:

- **SLURM** — HPC job scheduler that allocates compute resources and queues jobs on the cluster
- **Nextflow** — Workflow manager that orchestrates and parallelizes pipeline steps, handles container binding, and manages outputs
- **Singularity/Apptainer** — Container system that packages each bioinformatics tool with all its dependencies for reproducible execution

## How Nextflow Was Used

Nextflow orchestrates the variant calling workflow by running Clair3 and DeepVariant as independent parallel processes. Each tool is defined as a Nextflow `process` block in `variant_calling.nf` and executed inside its own Singularity container.

The workflow was submitted to SLURM using `nextflow_pipeline.sh`, which allocates 8 CPUs and 16GB RAM. Nextflow then manages the execution of each process, handles container binding via `nextflow.config`, and copies outputs to the results directory using `publishDir`.

Key Nextflow features used:
- **DSL2 syntax** — modern Nextflow workflow definition
- **process containers** — each tool runs in its own Singularity image
- **publishDir** — automatically copies outputs to the results folder
- **parallel execution** — Clair3 and DeepVariant run simultaneously
- **SLURM integration** — job submitted and managed by the HPC scheduler

To verify the SLURM + Nextflow + Singularity integration, the built-in Nextflow hello world pipeline was also successfully executed, confirming the full stack works correctly on this cluster.

## Pipeline Steps

1. **Read Alignment** — Raw PacBio HiFi reads aligned to GRCh38 using Minimap2 (map-hifi preset)
2. **BAM Processing** — Alignment sorted and indexed using Samtools
3. **Variant Calling with Clair3** — Neural network-based variant caller using hifi_revio model
4. **Variant Calling with DeepVariant** — Google's deep learning variant caller (PACBIO model)
5. **Chr1-22 Filtering** — VCF outputs filtered to chromosomes 1-22 using BCFtools
6. **Benchmarking** — Both VCFs evaluated against GIAB HG002 v4.2.1 truth set using hap.py

## How to Run

### Submit the full pipeline via SLURM + Nextflow + Singularity:
sbatch nextflow_pipeline.sh
### Or run the pipeline directly via SLURM:
sbatch pipeline.sh
### Monitor job status:
squeue -u your_username
## Benchmarking Results (Chr1-22)

Benchmarked against GIAB HG002 NISTv4.2.1 truth set using hap.py with vcfeval engine.
Note: A quarter subset of the full HG002 dataset was used. Low recall is expected due to limited read coverage across the genome. At full coverage, both tools typically achieve over 99% recall.

### SNP Performance

| Metric | Clair3 | DeepVariant |
|--------|--------|-------------|
| Variants Called | 252,391 | 163,639 |
| True Positives | 167,468 | 101,556 |
| False Positives | 35,508 | 31,719 |
| Precision | 82.5% | 76.2% |
| Recall | 4.98% | 3.02% |
| F1 Score | 0.094 | 0.058 |

### INDEL Performance

| Metric | Clair3 | DeepVariant |
|--------|--------|-------------|
| Variants Called | 26,778 | 24,264 |
| True Positives | 11,148 | 10,792 |
| False Positives | 5,164 | 5,825 |
| Precision | 68.8% | 64.9% |
| Recall | 2.12% | 2.05% |
| F1 Score | 0.041 | 0.040 |

### Key Findings

- Clair3 outperformed DeepVariant at this coverage level, achieving higher recall and precision for both SNPs and INDELs
- DeepVariant called fewer variants overall (163,639 vs 252,391 SNPs) reflecting its more conservative calling strategy
- Both tools maintained high precision despite low coverage, meaning calls made were mostly correct
- Low recall is attributed to the quarter-subset dataset — many genomic regions lacked sufficient read depth

## Repository Contents

| File | Description |
|------|-------------|
| variant_calling.nf | Nextflow workflow — runs Clair3 and DeepVariant in parallel |
| nextflow_pipeline.sh | SLURM job script that executes the Nextflow pipeline |
| pipeline.sh | SLURM job script for direct Apptainer execution |
| nextflow.config | Nextflow configuration for Singularity containers |
| deepvariant.nf | Standalone Nextflow script for DeepVariant |

## Tools and Containers

| Tool | Version | Container | Purpose |
|------|---------|-----------|---------|
| Minimap2 | latest | docker://staphb/minimap2 | Read alignment |
| Samtools | latest | docker://staphb/samtools | BAM processing |
| Clair3 | latest | docker://hkubal/clair3 | Variant calling |
| DeepVariant | 1.6.1 | docker://google/deepvariant:1.6.1 | Variant calling |
| BCFtools | latest | docker://staphb/bcftools | VCF filtering |
| hap.py | latest | docker://pkrusche/hap.py | Benchmarking |

## Data

- Sample: HG002 (Ashkenazi son, GIAB reference sample)
- Sequencing: PacBio HiFi (Revio chemistry)
- Input: 25% subset of full dataset (~502MB FASTQ)
- Reference: GRCh38 (hg38)
- Truth set: GIAB NISTv4.2.1 (chr1-22, GRCh38)
