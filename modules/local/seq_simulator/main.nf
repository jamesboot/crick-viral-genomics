process SEQ_SIMULATOR {
    tag "$meta.id"

    container "docker.io/thecrick/pipetech_seqsim:dev"

    input:
    tuple val(meta), path(refs)
    path(config)
    val(profile)
    val(num_reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('seq-sim'), eval('seq-sim --version'), topic: versions 

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seq-sim \\
    --hide-progress \\
    --log-file ${prefix}.seq-sim.log \\
    generate_reads \\
    -r $refs \\
    -c $config \\
    -p $profile \\
    -n $num_reads \\
    -o ${prefix}.fastq.gz
    """
}
