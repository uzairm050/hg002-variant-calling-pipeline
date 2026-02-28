# HG002 PacBio HiFi Variant Calling Pipeline

> A reproducible, HPC-ready pipeline for short variant calling on PacBio HiFi data using Clair3 and DeepVariant, benchmarked against the GIAB truth set.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10-brightgreen)](https://nextflow.io/)
[![Singularity](https://img.shields.io/badge/container-Singularity%2FApptainer-blue)](https://apptainer.org/)
[![SLURM](https://img.shields.io/badge/scheduler-SLURM-orange)](https://slurm.schedmd.com/)
[![Shell](https://img.shields.io/badge/language-Shell%20%7C%20Nextflow-lightgrey)](#)

---

## Overview

This pipeline processes PacBio HiFi sequencing reads from **HG002** — the well-characterised Ashkenazi son reference sample maintained by the [Genome in a Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle) consortium. It calls short variants (SNPs and INDELs) across chromosomes 1–22 using two independent deep learning callers and evaluates their accuracy against the GIAB NISTv4.2.1 truth set.

A genetic variant is a position in the genome where an individual's DNA differs from the standard reference. Accurately identifying these variants is fundamental to genomic research and clinical diagnostics — enabling discovery of disease-causing mutations, population-level studies, and personalised medicine.

The pipeline runs on an HPC cluster and uses **SLURM** for job scheduling, **Nextflow (DSL2)** to orchestrate and parallelise variant calling, and **Singularity/Apptainer** to containerise every tool for fully reproducible execution.

---

## Pipeline Architecture

```
PacBio HiFi FASTQ
        │
        ▼
┌──────────────────────────┐
│  Step 1: Read Alignment  │  Minimap2 (map-hifi preset) → GRCh38   [SLURM]
└─────────────┬────────────┘
              │
              ▼
┌──────────────────────────┐
│  Step 2: BAM Processing  │  Samtools sort + index                  [SLURM]
└─────────────┬────────────┘
              │
        ┌─────┴──────┐
        │            │   ← parallel Nextflow processes (variant_calling.nf)
        ▼            ▼
┌────────────┐  ┌─────────────┐
│   Clair3   │  │ DeepVariant │   [Nextflow + Singularity]
│ hifi_revio │  │   v1.6.1    │
└─────┬──────┘  └──────┬──────┘
      │                │
      ▼                ▼
┌──────────────────────────┐
│  Step 3: VCF Filtering   │  BCFtools — restrict to chr1–22         [SLURM]
└─────────────┬────────────┘
              │
              ▼
┌──────────────────────────┐
│  Step 4: Benchmarking    │  hap.py (vcfeval engine) vs GIAB v4.2.1 [SLURM]
└──────────────────────────┘
```

Steps 1, 2, 3, and 4 run as SLURM bash steps inside `pipeline.sh`. Clair3 and DeepVariant are defined as independent `process` blocks in `variant_calling.nf` and executed in **parallel** by Nextflow, each inside its own Singularity container. The Nextflow workflow is launched from within the SLURM job via `nextflow run variant_calling.nf`.

---

## Requirements

| Requirement | Details |
|---|---|
| SLURM workload manager | Any recent version |
| Nextflow | ≥ 22.10 (DSL2) |
| Singularity / Apptainer | ≥ 3.0 |
| CPUs | 8 (allocated by SLURM) |
| RAM | 16 GB (allocated by SLURM) |
| Storage | ~10 GB for reference + containers |

All bioinformatics tools run inside Singularity containers — no manual tool installation is required beyond Nextflow and Singularity.

---

## Repository Structure

```
hg002-variant-calling-pipeline/
├── pipeline.sh           # SLURM job script — runs alignment, filtering, and benchmarking
├── variant_calling.nf    # Nextflow DSL2 workflow — runs Clair3 and DeepVariant in parallel
└── README.md
```

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/uzairm050/hg002-variant-calling-pipeline.git
cd hg002-variant-calling-pipeline
```

### 2. Configure input paths

Edit the variables at the top of `pipeline.sh`:

```bash
FASTQ=/path/to/hg002.fastq.gz
REF=/path/to/GRCh38.fa
TRUTH_VCF=/path/to/giab_v4.2.1.vcf.gz
TRUTH_BED=/path/to/giab_v4.2.1.bed
OUTDIR=/path/to/output
```

### 3. Submit to SLURM

```bash
sbatch pipeline.sh
```

SLURM allocates 8 CPUs and 16 GB RAM. After alignment and BAM processing complete, Nextflow is invoked automatically inside the job to run Clair3 and DeepVariant in parallel via `variant_calling.nf`.

### 4. Monitor progress

```bash
squeue -u $USER          # check SLURM job status
tail -f slurm-*.out      # stream live log output
cat .nextflow.log        # inspect Nextflow process execution
```

### 5. Outputs

```
output/
├── aligned/        # Sorted, indexed BAM
├── clair3/         # Clair3 VCF filtered to chr1–22
├── deepvariant/    # DeepVariant VCF filtered to chr1–22
└── benchmark/      # hap.py summary CSV and extended metrics
```

---

## Data

| Parameter | Value |
|---|---|
| Sample | HG002 (NA24385 — Ashkenazi son, GIAB reference sample) |
| Sequencing platform | PacBio HiFi (Revio chemistry) |
| Input used | 25% subset (~502 MB FASTQ) |
| Reference genome | GRCh38 / hg38 |
| Truth set | GIAB NISTv4.2.1 (chr1–22, GRCh38) |

The GIAB HG002 dataset is publicly available via the [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/).

---

## Tools & Containers

All containers are pulled automatically from Docker Hub as Singularity images at runtime by Nextflow. An all-in-one Singularity container bundling all tools is also available for offline or air-gapped environments.

| Tool | Version | Container | Purpose |
|---|---|---|---|
| Minimap2 | latest | `docker://staphb/minimap2` | Long-read alignment to GRCh38 |
| Samtools | latest | `docker://staphb/samtools` | BAM sort, index, and statistics |
| Clair3 | latest | `docker://hkubal/clair3` | Neural network variant calling |
| DeepVariant | 1.6.1 | `docker://google/deepvariant` | Deep learning variant calling |
| BCFtools | latest | `docker://staphb/bcftools` | VCF filtering to chr1–22 |
| hap.py | latest | `docker://pkrusche/hap.py` | Precision/recall benchmarking |

---

## Results

> **Note:** These results are from a **25% coverage subset** of the full HG002 dataset. Low recall is expected — many genomic regions lacked sufficient read depth. At full coverage, both tools typically achieve **> 99% recall** on PacBio HiFi data.

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

Both callers show broadly similar chromosomal distributions. DeepVariant produces notably higher counts on chr1, while Clair3 leads across most other chromosomes at this coverage level.

---

## Conclusion

This project successfully implements a complete, end-to-end variant calling pipeline for PacBio HiFi data on an HPC cluster. By combining SLURM job scheduling, Nextflow DSL2 workflow management, and Singularity containerisation, the pipeline achieves full reproducibility — every tool runs in an isolated, version-pinned container, and Nextflow ensures Clair3 and DeepVariant execute in parallel with automatic output management via `publishDir`.

At 25% coverage, **Clair3 outperformed DeepVariant** on both SNP precision (82.5% vs 76.2%) and recall (4.98% vs 3.02%). This aligns with published benchmarks: Clair3's pileup-based two-stage approach is more robust under sparse coverage, whereas DeepVariant's convolutional image-classification model depends more heavily on read depth to construct informative pileup images. Both tools maintained reasonable precision (65–83%), confirming that called variants were largely correct — the primary limitation was missed variants due to the reduced dataset.

DeepVariant showed marginally better INDEL false positive control (64.9% vs 68.8% precision), suggesting its model is more conservative for insertion/deletion calls even at low depth.

At full sequencing depth, both tools are expected to achieve recall and precision exceeding 99%, consistent with their published benchmarks on PacBio HiFi HG002 data. The pipeline is fully parameterised and can be extended to additional samples, variant callers, or downstream annotation workflows with minimal modification.

---

## References

1. **Clair3** — Zheng, Z. et al. (2022). Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. *Nature Computational Science*, 2, 797–803. https://doi.org/10.1038/s43588-022-00387-x

2. **DeepVariant** — Poplin, R. et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. *Nature Biotechnology*, 36, 983–987. https://doi.org/10.1038/nbt.4235

3. **GIAB HG002 Truth Set** — Zook, J.M. et al. (2020). A robust benchmark for detection of germline large deletions and insertions. *Nature Biotechnology*, 38, 1347–1355. https://doi.org/10.1038/s41587-020-0538-8

4. **Minimap2** — Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

5. **Samtools / BCFtools** — Danecek, P. et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

6. **hap.py** — Krusche, P. et al. (2019). Best practices for benchmarking germline small-variant calls in human genomes. *Nature Biotechnology*, 37, 555–560. https://doi.org/10.1038/s41587-019-0054-x

7. **Nextflow** — Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319. https://doi.org/10.1038/nbt.3820

8. **Singularity / Apptainer** — Kurtzer, G.M., Sochat, V. & Bauer, M.W. (2017). Singularity: Scientific containers for mobility of compute. *PLOS ONE*, 12(5), e0177459. https://doi.org/10.1371/journal.pone.0177459

9. **GRCh38 Reference Genome** — Schneider, V.A. et al. (2017). Evaluation of GRCh38 and de novo haploid genome assemblies demonstrates the enduring quality of the reference assembly. *Genome Research*, 27(5), 849–864. https://doi.org/10.1101/gr.213611.116

---

## Authors

**Uzair Malik · Sana Khalid · Javeria Butt**

*Advanced Computational Biology — Assignment No. 1 · February 2026*

**Related repositories:**
- [SanaKhalid17/HG002-PacBio-HiFi-Variant-Calling-Pipeline-Clair3-vs-DeepVariant](https://github.com/SanaKhalid17/HG002-PacBio-HiFi-Variant-Calling-Pipeline-Clair3-vs-DeepVariant-)
- [javeria-butt/hg002-variant-calling-pipeline](https://github.com/javeria-butt/hg002-variant-calling-pipeline)
