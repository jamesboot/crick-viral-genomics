process SPADES_ASSEMBLE {
    tag "$meta.id"
    label "process_high"

    container "docker.io/thecrick/pipetech_spades:amd64-4.0.0"

    input:
    tuple val(meta), path(reads)

    output:
    // tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('spades'), eval('spades.py --version'), topic: versions 

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
    """
}
