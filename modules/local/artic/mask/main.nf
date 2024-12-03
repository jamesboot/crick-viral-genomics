process ARTIC_MASK {
    label 'process_single'
    tag "$meta.id"

    container "community.wave.seqera.io/library/artic_samtools:23997b5f9b2e3c39"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(mask)
    tuple val(meta3), path(fasta), path(fai)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    artic_mask \\
        $args \\
        $fasta \\
        $mask \\
        $vcf \\
        ${prefix}.preconsensus.fasta
    """
}
