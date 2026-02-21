#!/bin/bash
#SBATCH --job-name=variant_calling
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

WORKDIR=$(pwd)

echo "=== Step 1: Align reads ==="
apptainer exec --bind $WORKDIR:$WORKDIR $WORKDIR/minimap2_latest.sif \
  minimap2 -t 8 -ax map-hifi \
  $WORKDIR/data/GRCh38.fa \
  $WORKDIR/data/HG002_quarter.fastq.gz > $WORKDIR/results/aligned.sam

echo "=== Step 2: Sort BAM ==="
apptainer exec --bind $WORKDIR:$WORKDIR $WORKDIR/samtools_latest.sif \
  samtools sort -@ 8 -o $WORKDIR/results/aligned.bam $WORKDIR/results/aligned.sam

echo "=== Step 3: Index BAM ==="
apptainer exec --bind $WORKDIR:$WORKDIR $WORKDIR/samtools_latest.sif \
  samtools index $WORKDIR/results/aligned.bam

rm $WORKDIR/results/aligned.sam

echo "=== Step 4: Clair3 Variant Calling ==="
apptainer exec --bind $WORKDIR:$WORKDIR $WORKDIR/clair3_latest.sif \
  run_clair3.sh \
  --bam_fn=$WORKDIR/results/aligned.bam \
  --ref_fn=$WORKDIR/data/GRCh38.fa \
  --threads=8 --platform=hifi \
  --model_path=/opt/models/hifi_revio \
  --output=$WORKDIR/results/clair3_output

echo "=== Step 5: DeepVariant ==="
apptainer exec --bind $WORKDIR:$WORKDIR $WORKDIR/deepvariant_1.6.1.sif \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=$WORKDIR/data/GRCh38.fa \
  --reads=$WORKDIR/results/aligned.bam \
  --output_vcf=$WORKDIR/results/deepvariant_output.vcf.gz \
  --num_shards=8

echo "=== Step 6: Benchmark Clair3 ==="
apptainer exec --bind $WORKDIR:$WORKDIR $WORKDIR/hap.py_latest.sif \
  /opt/hap.py/bin/hap.py \
  $WORKDIR/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  $WORKDIR/results/clair3_output/merge_output.vcf.gz \
  -f $WORKDIR/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r $WORKDIR/data/GRCh38.fa \
  -o $WORKDIR/results/happy_clair3/happy \
  --engine=vcfeval --threads=8

echo "=== Step 7: Benchmark DeepVariant ==="
apptainer exec --bind $WORKDIR:$WORKDIR $WORKDIR/hap.py_latest.sif \
  /opt/hap.py/bin/hap.py \
  $WORKDIR/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  $WORKDIR/results/deepvariant_output.vcf.gz \
  -f $WORKDIR/data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r $WORKDIR/data/GRCh38.fa \
  -o $WORKDIR/results/happy_deepvariant/happy \
  --engine=vcfeval --threads=8

echo "=== Pipeline complete ==="
