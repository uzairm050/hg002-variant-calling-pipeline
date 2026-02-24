nextflow.enable.dsl=2

/*
Pipeline: FASTA + HiFi FASTQ -> minimap2 (map-hifi) -> samtools sort/index -> Clair3 -> VCF
Inputs:
  --reference     reference FASTA (GRCh38)
  --reads         HiFi reads FASTQ/FASTQ.GZ (quarter subset of HG002)
  --model_path    Clair3 model directory for HiFi (hifi_revio)
Outputs:
  outdir/<sample>.sorted.bam(.bai) + outdir/<sample>.vcf.gz(.tbi)
*/

if( !params.reference )  error "Missing --reference"
if( !params.reads )      error "Missing --reads"
if( !params.model_path ) error "Missing --model_path"

Channel.fromPath(params.reference).set { ch_ref }
Channel.fromPath(params.reads).set { ch_reads }

def inferSampleName(readsFile) {
  def base = readsFile.getBaseName()
  base = base.replaceAll(/(\.fastq|\.fq)$/,'')
  return base
}

def minimap2_container = params.use_local_sifs ? params.sif_minimap2 : "docker://quay.io/biocontainers/minimap2:2.28--he4a0461_0"
def samtools_container = params.use_local_sifs ? params.sif_samtools : "docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
def clair3_container   = params.use_local_sifs ? params.sif_clair3   : "docker://hkubal/clair3:latest"

workflow {
  def sampleName = params.sample ?: inferSampleName(file(params.reads))
  SETUP_CHECK(ch_ref, ch_reads, sampleName)
  def ch_mmi = INDEX_REF_MINIMAP2(ch_ref)
  def ch_sam = ALIGN_MINIMAP2(ch_reads, ch_mmi, ch_ref, sampleName)
  def ch_bam = SORT_INDEX_SAMTOOLS(ch_sam, sampleName)
  CALL_VARIANTS_CLAIR3(ch_bam, ch_ref, sampleName)
}

process SETUP_CHECK {
  tag { sample }
  container minimap2_container
  input:
    path ref_fa
    path reads
    val  sample
  output:
    val(true)
  script:
    """
    set -euo pipefail
    echo "=== SETUP CHECK ==="
    echo "Sample: ${sample}"
    echo "Reference: ${ref_fa}"
    echo "Reads: ${reads}"
    echo "Platform: ${params.platform}"
    echo "Minimap2 preset: ${params.minimap2_preset}"
    echo ""
    test -s "${ref_fa}" || { echo "ERROR: reference FASTA is missing/empty"; exit 1; }
    test -s "${reads}"  || { echo "ERROR: reads FASTQ is missing/empty"; exit 1; }
    if [ ! -d "${params.model_path}" ]; then
      echo "ERROR: Clair3 model_path not found: ${params.model_path}" >&2
      exit 1
    fi
    echo "Tool versions:"
    minimap2 --version | head -n 1 || true
    echo "Input checks passed."
    """
}

process INDEX_REF_MINIMAP2 {
  tag "ref_index"
  container minimap2_container
  input:
    path ref_fa
  output:
    path "reference.mmi", emit: mmi
  script:
    """
    set -euo pipefail
    minimap2 -d reference.mmi ${ref_fa}
    """
}

process ALIGN_MINIMAP2 {
  tag { sample }
  publishDir params.outdir, mode: 'copy'
  container minimap2_container
  input:
    path reads
    path mmi
    path ref_fa
    val  sample
  output:
    path "${sample}.sam", emit: sam
  script:
    """
    set -euo pipefail
    minimap2 -t ${task.cpus} -a -x ${params.minimap2_preset} ${mmi} ${reads} > ${sample}.sam
    """
}

process SORT_INDEX_SAMTOOLS {
  tag { sample }
  publishDir params.outdir, mode: 'copy'
  container samtools_container
  input:
    path sam
    val  sample
  output:
    path "${sample}.sorted.bam", emit: bam
    path "${sample}.sorted.bam.bai", emit: bai
  script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${sam}
    samtools index -@ ${task.cpus} ${sample}.sorted.bam
    """
}

process CALL_VARIANTS_CLAIR3 {
  tag { sample }
  publishDir params.outdir, mode: 'copy'
  container clair3_container
  input:
    path bam
    path bai
    path ref_fa
    val  sample
  output:
    path "${sample}.vcf.gz"
    path "${sample}.vcf.gz.tbi"
  script:
    """
    set -euo pipefail
    mkdir -p clair3_out
    run_clair3.sh \
      --bam_fn ${bam} \
      --ref_fn ${ref_fa} \
      --threads ${task.cpus} \
      --platform ${params.platform} \
      --model_path ${params.model_path} \
      --output clair3_out \
      ${params.clair3_args}
    if [ -f clair3_out/merge_output.vcf.gz ]; then
      cp clair3_out/merge_output.vcf.gz ${sample}.vcf.gz
      cp clair3_out/merge_output.vcf.gz.tbi ${sample}.vcf.gz.tbi
    else
      vcf=\$(find clair3_out -maxdepth 4 -name "*.vcf.gz" | head -n 1)
      if [ -z "\$vcf" ]; then
        echo "ERROR: Clair3 did not produce a VCF." >&2
        exit 1
      fi
      cp "\$vcf" ${sample}.vcf.gz
      if [ -f "\$vcf.tbi" ]; then
        cp "\$vcf.tbi" ${sample}.vcf.gz.tbi
      elif command -v tabix >/dev/null 2>&1; then
        tabix -p vcf ${sample}.vcf.gz
      else
        echo "WARNING: No VCF index produced." >&2
      fi
    fi
    """
}
