process MEDAKA_ANNOTATE {
    label 'process_single'
    tag "$meta.id"

    // medaka v2.0.0
    container "docker.io/ontresearch/medaka:sha447c70a639b8bcf17dc49b51e74dfcde6474837b"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(bam), path(bai)
    val read_group

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rg = read_group ? "--RG ${read_group[0]}" : ""
    """
    medaka tools annotate \\
        $args \\
        $rg \\
        $vcf \\
        $bam \\
        $fasta \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
