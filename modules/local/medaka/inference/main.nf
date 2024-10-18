process MEDAKA_INFERENCE {
    label 'process_medium'
    tag "$meta.id"

    // medaka v2.0.0
    container "docker.io/ontresearch/medaka:sha447c70a639b8bcf17dc49b51e74dfcde6474837b"

    input:
    tuple val(meta), path(bam), path(bai)
    val read_group

    output:
    tuple val(meta), path("*.hdf"), emit: hdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rg = read_group ? "--RG ${read_group[0]}" : ""
    """
    medaka inference \\
        $args \\
        --threads $task.cpus \\
        $rg \\
        $bam \\
        ${prefix}.hdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
