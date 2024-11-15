process GFF_FLU {
    tag "$meta.id"
    label "process_single"

    container "community.wave.seqera.io/library/gfflu:0.0.2--bbd82e111d01609a"

    input:
    tuple val(meta), path(fasta)

    // output:
    // tuple val(meta), path("*.vcf"), emit: vcf
    // tuple val("${task.process}")  , val('LoFreq'), eval('lofreq version'), topic: versions

    script:
    // def args   = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"

    """
    lofreq call-parallel \\
        --pp-threads $task.cpus \\
        -f $fasta \\
        $args \\
        -o ${prefix}.vcf \\
        $bam
    """
}
