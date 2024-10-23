process ARTIC_MAKE_DEPTH_MASK {
    label 'process_single'
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/artic_samtools:add476d24521086e':
        'community.wave.seqera.io/library/artic_samtools:23997b5f9b2e3c39' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path('*.txt')   , emit: mask
    tuple val(meta), path('*.depths'), emit: depths, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    artic_make_depth_mask \\
        $args \\
        $fasta \\
        $bam \\
        ${prefix}.coverage_mask.txt
    """
}
