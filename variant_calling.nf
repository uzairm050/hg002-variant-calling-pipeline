#!/usr/bin/env nextflow

params.bam    = "${PWD}/results/aligned.bam"
params.ref    = "${PWD}/data/GRCh38.fa"
params.outdir = "${PWD}/results"
params.model  = "/opt/models/hifi_revio"

process CLAIR3 {
    container "${PWD}/clair3_latest.sif"
    publishDir "${params.outdir}/clair3_nf", mode: 'copy'

    output:
    path "clair3_out/*"

    script:
    """
    mkdir -p clair3_out
    run_clair3.sh \
        --bam_fn=${params.bam} \
        --ref_fn=${params.ref} \
        --threads=8 \
        --platform=hifi \
        --model_path=${params.model} \
        --output=\${PWD}/clair3_out
    """
}

process DEEPVARIANT {
    container "${PWD}/deepvariant_1.6.1.sif"
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path "deepvariant_nf.vcf.gz"

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref=${params.ref} \
        --reads=${params.bam} \
        --output_vcf=\${PWD}/deepvariant_nf.vcf.gz \
        --num_shards=8
    """
}

workflow {
    CLAIR3()
    DEEPVARIANT()
}
