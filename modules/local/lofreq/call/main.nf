process LOFREQ_CALL {
    tag "$meta.id"
    label "process_medium"

    container "docker.io/thecrick/pipetech_lowfreq:amd64-2.1.5"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val("${task.process}")  , val('LoFreq'), eval('lofreq version'), topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    lofreq call-parallel \\
        --pp-threads $task.cpus \\
        -f $fasta \\
        $args \\
        -o ${prefix}.vcf \\
        $bam
    """
}
