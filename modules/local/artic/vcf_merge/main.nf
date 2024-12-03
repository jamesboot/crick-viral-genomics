process ARTIC_VCF_MERGE {
    label 'process_single'
    tag "$meta.id"

    container "community.wave.seqera.io/library/artic_samtools:23997b5f9b2e3c39"

    input:
    tuple val(meta), path(vcf), val(pool)
    path bed

    output:
    tuple val(meta), path('*.merged.vcf') , emit: vcf
    tuple val(meta), path('*.primers.vcf'), emit: primer_vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_pool = ""
    vcf.eachWithIndex { item, index ->
        vcf_pool += " ${pool[index]}:${item}"
    }
    """
    artic_vcf_merge \\
        $prefix \\
        $bed \\
        $vcf_pool
    """
}
