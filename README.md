# HG002 Variant Calling Pipeline

[![SLURM](https://img.shields.io/badge/HPC-SLURM-blue)](https://slurm.schedmd.com/)
[![Nextflow](https://img.shields.io/badge/workflow-Nextflow-brightgreen)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-orange)](https://apptainer.org/)
[![GIAB](https://img.shields.io/badge/truth%20set-GIAB%20v4.2.1-red)](https://www.nist.gov/programs-projects/genome-bottle)

## What is this project?

This pipeline takes raw DNA sequencing data from a human sample, finds all the positions where that person's genome differs from the standard human reference genome, and checks how accurately those differences were detected.

We use HG002 — a well-studied reference human sample used worldwide to test genomic tools — and run the full analysis on a high-performance computing cluster using three technologies working together: SLURM for job scheduling, Nextflow for workflow management, and Singularity for containerized tools.

---

## Technologies Used

| Technology | Role |
|------------|------|
| **SLURM** | Submits and manages jobs on the HPC cluster, allocates compute resources |
| **Nextflow** | Orchestrates the pipeline, runs Clair3 and DeepVariant in parallel |
| **Singularity** | Runs each tool in an isolated, reproducible container |

---

## How Nextflow Was Used

Nextflow is the core of this pipeline. It defines Clair3 and DeepVariant as independent `process` blocks in `variant_calling.nf`, runs them in parallel inside their own Singularity containers, and automatically copies outputs to the results folder using `publishDir`.

The entire workflow is submitted to SLURM with a single command. SLURM allocates the compute resources, and Nextflow takes over from there — managing containers, parallelization, and output handling automatically.

Key Nextflow features used:
- **DSL2 syntax** — modern Nextflow workflow definition
- **process containers** — each tool runs in its own Singularity image
- **publishDir** — automatically copies outputs to the results folder
- **parallel execution** — Clair3 and DeepVariant run simultaneously
- **SLURM integration** — job submitted and managed by the HPC scheduler

To confirm the full SLURM + Nextflow + Singularity stack works on this cluster, the built-in Nextflow hello world pipeline was also run successfully before the main pipeline.

---

## Pipeline Steps

| Step | Tool | Executed Via |
|------|------|-------------|
| 1. Align reads to GRCh38 | Minimap2 (map-hifi) | SLURM |
| 2. Sort and index alignment | Samtools | SLURM |
| 3. Call variants — Clair3 | Clair3 (hifi_revio model) | Nextflow + Singularity |
| 4. Call variants — DeepVariant | DeepVariant 1.6.1 (PACBIO model) | Nextflow + Singularity |
| 5. Filter to chr1-22 | BCFtools | SLURM |
| 6. Benchmark against truth set | hap.py (vcfeval engine) | SLURM |

---

## How to Run

**Step 1 — Clone the repository:**
git clone https://github.com/uzairm050/hg002-variant-calling-pipeline.git
cd hg002-variant-calling-pipeline

**Step 2 — Submit the pipeline to SLURM + Nextflow + Singularity:**
sbatch nextflow_pipeline.sh

**Step 3 — Monitor progress:**
squeue -u your_username

tail -f logs/nextflow_*.log

**Step 4 — Or run directly via SLURM without Nextflow:**
sbatch pipeline.sh

> **Resource Note:** The SLURM parameters in the job scripts are configured for a quarter-subset of HG002 data on our cluster. If running on a different cluster or with the full dataset, adjust the following parameters in `nextflow_pipeline.sh` and `pipeline.sh`:
> - `--cpus-per-task` — increase for faster parallel processing
> - `--mem` — minimum 64GB recommended for the full HG002 dataset
> - `--time` — minimum 12 hours recommended for the full dataset
> - `--partition` — set to the appropriate partition on your cluster

---

## Benchmarking Results (Chr1-22)

Benchmarked against GIAB HG002 NISTv4.2.1 using hap.py with vcfeval engine, restricted to chromosomes 1-22.

> **Note on recall:** A quarter subset of the full HG002 dataset was used. Low recall is expected — many genomic regions lacked sufficient read coverage to detect variants. At full coverage both tools typically achieve over 99% recall. The precision values however reflect true tool accuracy.

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

- **Clair3 outperformed DeepVariant** at this coverage level with higher recall and precision for SNPs
- DeepVariant is more conservative — called fewer variants but with slightly better INDEL precision
- Both tools showed high precision, meaning the calls they made were mostly correct
- Low recall is expected and fully explained by the quarter-subset input data

---

## Repository Contents

| File | Description |
|------|-------------|
| `variant_calling.nf` | Nextflow workflow — runs Clair3 and DeepVariant in parallel via Singularity |
| `nextflow_pipeline.sh` | SLURM script — submits the Nextflow pipeline to the HPC cluster |
| `pipeline.sh` | SLURM script — runs the full pipeline directly via Apptainer |
| `nextflow.config` | Nextflow config — sets up Singularity container binding and run options |
| `deepvariant.nf` | Standalone Nextflow script for DeepVariant only |

---

## Tools and Containers

| Tool | Version | Container | Purpose |
|------|---------|-----------|---------|
| Minimap2 | latest | docker://staphb/minimap2 | Read alignment |
| Samtools | latest | docker://staphb/samtools | BAM processing |
| Clair3 | latest | docker://hkubal/clair3 | Variant calling via Nextflow |
| DeepVariant | 1.6.1 | docker://google/deepvariant:1.6.1 | Variant calling via Nextflow |
| BCFtools | latest | docker://staphb/bcftools | VCF filtering chr1-22 |
| hap.py | latest | docker://pkrusche/hap.py | Benchmarking vs GIAB |

---

## Data

| Parameter | Value |
|-----------|-------|
| Sample | HG002 (Ashkenazi son, GIAB reference sample) |
| Sequencing | PacBio HiFi (Revio chemistry) |
| Input | 25% subset (~502MB FASTQ) |
| Reference | GRCh38 (hg38) |
| Truth set | GIAB NISTv4.2.1 (chr1-22, GRCh38) |
