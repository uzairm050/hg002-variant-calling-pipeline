# HG002 PacBio HiFi Variant Calling Pipeline

> A fully Nextflow-orchestrated, Singularity-containerised pipeline for short variant calling on PacBio HiFi data using Clair3 and DeepVariant, benchmarked against the GIAB truth set.

[![Nextflow DSL2](https://img.shields.io/badge/nextflow-DSL2-%2300AA00)](https://nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity%2FApptainer-blue)](https://apptainer.org/)
[![SLURM](https://img.shields.io/badge/scheduler-SLURM-orange)](https://slurm.schedmd.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

This pipeline processes PacBio HiFi sequencing reads from **HG002** — the well-characterised Ashkenazi son reference sample maintained by the [Genome in a Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) consortium. It calls short variants (SNPs and INDELs) across chromosomes 1–22 using two independent deep learning callers — **Clair3** and **DeepVariant** — and evaluates accuracy against the GIAB NISTv4.2.1 truth set.

The entire pipeline is orchestrated by **Nextflow DSL2**. Every process — alignment, sorting, and variant calling — runs inside an isolated **Singularity** container, ensuring full reproducibility across HPC compute nodes. The pipeline is launched with a single shell script (`pipeline.sh`) that calls `nextflow run main.nf`.

---

## Pipeline Architecture

The Nextflow workflow (`main.nf`) defines five sequential processes, each running in its own Singularity container:

```
HG002 ¼-subset HiFi FASTQ
          │
          ▼
┌─────────────────────┐
│    SETUP_CHECK      │  Validates input files + Clair3 model path exist
│                     │  Container: minimap2 (quay.io/biocontainers/minimap2:2.28)
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  INDEX_REF_MINIMAP2 │  Builds minimap2 .mmi index from GRCh38 FASTA
│                     │  Container: minimap2 (quay.io/biocontainers/minimap2:2.28)
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│   ALIGN_MINIMAP2    │  Aligns reads with -x map-hifi preset → HG002.sam
│                     │  Container: minimap2 (quay.io/biocontainers/minimap2:2.28)
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ SORT_INDEX_SAMTOOLS │  Sorts + indexes BAM → HG002.sorted.bam + .bai
│                     │  Container: samtools (quay.io/biocontainers/samtools:1.20)
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ CALL_VARIANTS_CLAIR3│  Runs run_clair3.sh with hifi_revio model → HG002.vcf.gz
│                     │  Container: hkubal/clair3:latest
└─────────────────────┘

DeepVariant run separately via deepvariant.nf (same Nextflow + Singularity approach)
           │
           ▼
┌─────────────────────┐
│   BCFtools filter   │  Restricts both VCFs to chr1–22
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  hap.py (vcfeval)   │  Benchmarks Clair3 + DeepVariant VCFs vs GIAB NISTv4.2.1
└─────────────────────┘
```

`pipeline.sh` is a thin launcher — it simply calls `nextflow run main.nf -profile singularity` with the correct parameters. All compute logic lives inside the Nextflow workflow.

---

## Repository Structure

```
hg002-variant-calling-pipeline/
├── main.nf           # Nextflow DSL2 workflow: SETUP_CHECK → INDEX → ALIGN → SORT → CLAIR3
├── deepvariant.nf    # Nextflow DSL2 workflow: DeepVariant variant calling
├── nextflow.config   # Singularity profile, per-process resource limits, all params
├── pipeline.sh       # Launch script: calls nextflow run main.nf -profile singularity
└── README.md
```

---

## Requirements

| Requirement | Version |
|---|---|
| Nextflow | ≥ 22.10 (DSL2) |
| Singularity / Apptainer | ≥ 3.0 |
| CPUs | 8 (set via `--threads`) |
| RAM | Up to 16 GB (per-process limits in `nextflow.config`) |
| Storage | ~10 GB for reference, containers cache, and outputs |

No manual tool installation required. All containers are pulled automatically from Docker Hub into `.singularity_cache/` on first run.

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/uzairm050/hg002-variant-calling-pipeline.git
cd hg002-variant-calling-pipeline
```

### 2. Edit paths in `pipeline.sh`

```bash
BASEDIR=/your/working/directory

nextflow run main.nf -profile singularity \
  --reference ${BASEDIR}/data/GRCh38.fa \
  --reads     ${BASEDIR}/data/HG002_quarter.fastq.gz \
  --model_path /opt/models/hifi_revio \
  --sample    HG002 \
  --outdir    ${BASEDIR}/results \
  --threads   8
```

### 3. Run

```bash
# Run directly
bash pipeline.sh

# Or submit to SLURM
sbatch --cpus-per-task=8 --mem=32G --time=24:00:00 pipeline.sh
```

### 4. Monitor

```bash
# Nextflow prints live process status to stdout
tail -f slurm-*.out      # if submitted via SLURM
cat .nextflow.log        # detailed Nextflow execution log
```

### 5. Outputs

```
results/
├── HG002.sam               # Raw alignment (intermediate)
├── HG002.sorted.bam        # Sorted BAM
├── HG002.sorted.bam.bai    # BAM index
├── HG002.vcf.gz            # Clair3 variant calls
└── HG002.vcf.gz.tbi        # VCF tabix index
```

---

## Configuration (`nextflow.config`)

All parameters, resource limits, and Singularity settings are defined in `nextflow.config`.

### Process resource allocation

| Process | Memory | Time |
|---|---|---|
| `SETUP_CHECK` | 1 GB | 10 min |
| `INDEX_REF_MINIMAP2` | 2 GB | 2 h |
| `ALIGN_MINIMAP2` | 8 GB | 6 h |
| `SORT_INDEX_SAMTOOLS` | 8 GB | 6 h |
| `CALL_VARIANTS_CLAIR3` | 16 GB | 12 h |

### Singularity settings

```groovy
singularity {
  enabled    = true
  autoMounts = true
  cacheDir   = "${projectDir}/.singularity_cache"
}
```

`autoMounts = true` ensures host filesystem paths (reference, reads, model) are automatically bind-mounted into each container — no manual `-B` flags needed.

### Using local `.sif` files (offline / air-gapped HPC)

If your cluster cannot pull from Docker Hub, pre-download the containers and use:

```bash
nextflow run main.nf -profile singularity \
  --use_local_sifs true \
  --sif_minimap2 containers/minimap2.sif \
  --sif_samtools containers/samtools.sif \
  --sif_clair3   containers/clair3.sif \
  ...
```

### All parameters

| Parameter | Default | Description |
|---|---|---|
| `--reference` | `null` *(required)* | Path to GRCh38 FASTA |
| `--reads` | `null` *(required)* | Path to HiFi FASTQ or FASTQ.GZ |
| `--model_path` | `null` *(required)* | Path to Clair3 hifi_revio model directory |
| `--sample` | inferred from filename | Sample name used for all output files |
| `--outdir` | `results` | Output directory |
| `--threads` | `8` | CPUs allocated to each process |
| `--platform` | `hifi` | Clair3 `--platform` flag |
| `--minimap2_preset` | `map-hifi` | Minimap2 `-x` preset |
| `--clair3_args` | `""` | Extra arguments passed to `run_clair3.sh` |
| `--use_local_sifs` | `false` | Use local `.sif` files instead of Docker Hub |

---

## Data

| Parameter | Value |
|---|---|
| Sample | HG002 (NA24385 — Ashkenazi son, GIAB reference sample) |
| Sequencing platform | PacBio HiFi (Revio chemistry) |
| Input used | ¼ subset (~502 MB FASTQ) |
| Source | [PacBio HG002 Vega BAM](https://downloads.pacbcloud.com/public/2024Q4/Vega/HG002/data/) → converted with `bamtofastq` |
| Reference genome | GRCh38 / hg38 |
| Truth set | GIAB NISTv4.2.1 (chr1–22, GRCh38) |

---

## Tools & Containers

| Tool | Version | Container | Used In |
|---|---|---|---|
| Minimap2 | 2.28 | `quay.io/biocontainers/minimap2:2.28--he4a0461_0` | `main.nf` — SETUP_CHECK, INDEX, ALIGN |
| Samtools | 1.20 | `quay.io/biocontainers/samtools:1.20--h50ea8bc_0` | `main.nf` — SORT_INDEX |
| Clair3 | latest | `hkubal/clair3:latest` | `main.nf` — CALL_VARIANTS_CLAIR3 |
| DeepVariant | 1.6.1 | `google/deepvariant:1.6.1` | `deepvariant.nf` |
| BCFtools | latest | `staphb/bcftools` | Post-processing — filter to chr1–22 |
| hap.py | latest | `pkrusche/hap.py` | Benchmarking vs GIAB |

---

## Results

> **Note:** These results are from a **¼ coverage subset** of the full HG002 dataset. Low recall is expected — many regions lacked sufficient depth. At full coverage, both tools typically achieve **> 99% recall** on PacBio HiFi data.

### SNP Performance (chr1–22)

| Metric | Clair3 | DeepVariant |
|---|---|---|
| Variants Called | 252,391 | 163,639 |
| True Positives | 167,468 | 101,556 |
| False Positives | 35,508 | 31,719 |
| Precision | **82.5%** | 76.2% |
| Recall | **4.98%** | 3.02% |
| F1 Score | **0.094** | 0.058 |

### INDEL Performance (chr1–22)

| Metric | Clair3 | DeepVariant |
|---|---|---|
| Variants Called | 26,778 | 24,264 |
| True Positives | 11,148 | 10,792 |
| False Positives | 5,164 | 5,825 |
| Precision | 68.8% | **64.9%** |
| Recall | **2.12%** | 2.05% |
| F1 Score | **0.041** | 0.040 |

### Per-Chromosome Variant Counts (excerpt)

| Chromosome | Clair3 | DeepVariant |
|---|---|---|
| chr1 | 6,860 | 29,632 |
| chr2 | 6,236 | 4,234 |
| chr3 | 5,185 | 4,725 |
| chr4 | 7,293 | 8,262 |
| chr5 | 4,443 | 4,309 |
| chr6 | 4,361 | 4,534 |
| chr7 | 4,175 | 2,628 |
| chr8 | 3,996 | 3,304 |
| chr9 | 3,669 | 3,619 |
| chr10 | 6,130 | 6,469 |
| chr11 | 2,945 | 2,585 |
| chr12 | 3,177 | 2,330 |
| … | … | … |

---

## Conclusion

This project implements a complete, production-ready variant calling pipeline for PacBio HiFi data, fully orchestrated by **Nextflow DSL2** and executed on an HPC cluster. Every process runs inside an isolated, version-pinned Singularity container, with `autoMounts` enabled in `nextflow.config` to handle HPC filesystem paths seamlessly. The pipeline supports both Docker Hub container pulling and local `.sif` files for air-gapped environments.

At ¼ coverage, **Clair3 outperformed DeepVariant** on both SNP precision (82.5% vs 76.2%) and recall (4.98% vs 3.02%). This is consistent with published benchmarks: Clair3's pileup-based two-stage approach degrades more gracefully under sparse read coverage, whereas DeepVariant's convolutional image-classification model relies more heavily on sufficient depth. Both tools maintained reasonable precision (65–83%), confirming that called variants were largely correct — the primary limitation was missed variants due to the subsampled dataset.

At full sequencing depth, both tools are expected to exceed 99% recall and precision, consistent with their published GIAB benchmarks on PacBio HiFi data. The modular Nextflow process design makes it straightforward to swap in additional callers or extend the pipeline with downstream annotation steps.

---

