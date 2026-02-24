# HG002 Variant Calling Pipeline

[![Nextflow](https://img.shields.io/badge/workflow-Nextflow-brightgreen)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity-orange)](https://apptainer.org/)
[![GIAB](https://img.shields.io/badge/truth%20set-GIAB%20v4.2.1-red)](https://www.nist.gov/programs-projects/genome-bottle)

**Workflow:** minimap2 (map-hifi) → samtools sort/index → Clair3 → VCF → hap.py benchmarking

---

## What is this project?

This pipeline processes PacBio HiFi sequencing reads from HG002 — a well-studied reference human sample maintained by GIAB — aligns them to the GRCh38 reference genome, calls genetic variants using Clair3, and benchmarks the results against the GIAB v4.2.1 truth set.

---

## Learning Objectives

1. Run a reproducible genomics workflow using Nextflow and Singularity.
2. Align PacBio HiFi reads to a reference genome using minimap2.
3. Create a coordinate-sorted BAM and BAM index using samtools.
4. Call small variants (SNVs/indels) using Clair3.
5. Benchmark variant calls against a GIAB truth set using hap.py.

---

## Prerequisites

- **Nextflow** (DSL2)
- **Singularity** or **Apptainer**
- **Java** (required by Nextflow)

---

## Project Structure
```
hg002-variant-calling-pipeline/
├── main.nf               # Nextflow pipeline
├── nextflow.config       # Singularity and resource configuration
├── run_pipeline.sh       # Convenience wrapper script
├── README.md             # This file
├── data/
│   ├── GRCh38.fa                                          # reference genome
│   ├── HG002_quarter.fastq.gz                             # HiFi reads
│   ├── HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz         # GIAB truth set
│   └── HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
└── containers/
    ├── minimap2.sif
    ├── samtools.sif
    └── clair3.sif
```

---

## Inputs

### Required

| Parameter | Description |
|-----------|-------------|
| `--reference` | Reference FASTA (GRCh38) |
| `--reads` | HiFi FASTQ/FASTQ.GZ |
| `--model_path` | Clair3 HiFi model directory |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sample` | derived from FASTQ name | Output file label |
| `--outdir` | `results` | Output folder |
| `--threads` | `8` | CPU threads |

---

## How to Run

**Using the wrapper script:**
```bash
bash run_pipeline.sh
```

**Or directly with Nextflow:**
```bash
nextflow run main.nf -profile singularity \
  --reference data/GRCh38.fa \
  --reads data/HG002_quarter.fastq.gz \
  --model_path /shared/clair3_models/hifi \
  --sample HG002 \
  --outdir results \
  --threads 8
```

**Using local SIF files (recommended on HPC):**
```bash
nextflow run main.nf -profile singularity \
  --use_local_sifs true \
  --sif_minimap2 containers/minimap2.sif \
  --sif_samtools containers/samtools.sif \
  --sif_clair3   containers/clair3.sif \
  --reference data/GRCh38.fa \
  --reads data/HG002_quarter.fastq.gz \
  --model_path /shared/clair3_models/hifi \
  --sample HG002 \
  --outdir results \
  --threads 8
```

**Monitor progress:**
```bash
tail -f .nextflow.log
```

---

## Pipeline Steps

| Step | Process | Tool | Container |
|------|---------|------|-----------|
| 1 | `SETUP_CHECK` | Validates inputs and model path | minimap2 |
| 2 | `INDEX_REF_MINIMAP2` | Indexes reference genome | minimap2 |
| 3 | `ALIGN_MINIMAP2` | Aligns HiFi reads (map-hifi preset) | minimap2 |
| 4 | `SORT_INDEX_SAMTOOLS` | Sorts and indexes BAM | samtools |
| 5 | `CALL_VARIANTS_CLAIR3` | Calls SNVs/INDELs | clair3 |

---

## Expected Outputs
```
results/
├── HG002.sam
├── HG002.sorted.bam
├── HG002.sorted.bam.bai
├── HG002.vcf.gz
└── HG002.vcf.gz.tbi
```

---

## Benchmarking Results (Chr1-22)

Variants were filtered to chromosomes 1-22 and benchmarked against GIAB HG002 NISTv4.2.1 using hap.py with vcfeval engine.
```bash
hap.py \
  data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  results/HG002.vcf.gz \
  -f data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r data/GRCh38.fa \
  -o results/happy/happy \
  --engine=vcfeval \
  --threads=8
```

### SNP Performance

| Metric | Clair3 |
|--------|--------|
| Variants Called | 252,391 |
| True Positives | 167,468 |
| False Positives | 35,508 |
| Precision | 82.5% |
| Recall | 4.98% |
| F1 Score | 0.094 |

### INDEL Performance

| Metric | Clair3 |
|--------|--------|
| Variants Called | 26,778 |
| True Positives | 11,148 |
| False Positives | 5,164 |
| Precision | 68.8% |
| Recall | 2.12% |
| F1 Score | 0.041 |

> **Note:** A quarter subset of HG002 was used so low recall is expected. At full coverage Clair3 typically achieves over 99% recall on PacBio HiFi data.

---

## Tools and Containers

| Tool | Version | Container | Purpose |
|------|---------|-----------|---------|
| Minimap2 | 2.28 | quay.io/biocontainers/minimap2 | Read alignment |
| Samtools | 1.20 | quay.io/biocontainers/samtools | BAM processing |
| Clair3 | latest | docker://hkubal/clair3 | Variant calling |
| hap.py | latest | docker://pkrusche/hap.py | Benchmarking |

---

## Data

| Parameter | Value |
|-----------|-------|
| Sample | HG002 (Ashkenazi son, GIAB) |
| Sequencing | PacBio HiFi (Revio chemistry) |
| Input | 25% subset (~502MB FASTQ) |
| Reference | GRCh38 (hg38) |
| Truth set | GIAB NISTv4.2.1 (chr1-22, GRCh38) |
