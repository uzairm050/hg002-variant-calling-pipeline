# HG002 PacBio HiFi Variant Calling Pipeline

A variant calling pipeline for PacBio HiFi reads using Clair3 and DeepVariant, benchmarked against the GIAB truth set.

## Usage
Run with SLURM: sbatch pipeline.sh

## Pipeline Steps
1. Alignment with Minimap2 (map-hifi preset)
2. BAM sorting and indexing with Samtools
3. Variant calling with Clair3 (hifi_revio model)
4. Variant calling with DeepVariant 1.6.1 (PACBIO model)
5. Benchmarking with hap.py vs GIAB v4.2.1

## Results
| Metric | Clair3 | DeepVariant |
|--------|--------|-------------|
| Total Variants | 294,341 | 215,091 |
| SNP Recall | 4.98% | 3.02% |
| SNP Precision | 82.5% | 76.2% |
| INDEL Recall | 2.73% | 2.05% |
| INDEL Precision | 62.6% | 64.9% |

Note: Low recall expected - quarter subset of reads used.

## Containers
- Minimap2: docker://staphb/minimap2
- Samtools: docker://staphb/samtools
- Clair3: docker://hkubal/clair3
- DeepVariant: docker://google/deepvariant:1.6.1
- hap.py: docker://pkrusche/hap.py
