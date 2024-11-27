process SPADES_ASSEMBLE {
    tag "$meta.id"
    label "process_high"

    container "docker.io/thecrick/pipetech_spades:amd64-4.0.0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("spades/contigs.fasta"), emit: contigs
    path "versions.yml"                          , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    spades.py \\
        --threads $task.cpus \\
        --pe1-1 ${reads[0]} \\
        --pe1-2 ${reads[1]} \\
        $args \\
        -o ./spades

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$( spades.py --version 2>&1 | sed -n 's/^SPAdes genome assembler v//p' )
    END_VERSIONS
    """
}
