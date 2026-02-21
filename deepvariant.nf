process DEEPVARIANT {
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=PACBIO \\
        --ref=${params.ref} \\
        --reads=${params.bam} \\
        --output_vcf=${params.outdir}/deepvariant_output.vcf.gz \\
        --num_shards=8
    """
}

workflow {
    DEEPVARIANT()
}
