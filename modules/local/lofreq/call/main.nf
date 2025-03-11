process LOFREQ_CALL {
    tag "$meta.id"
    label "process_medium"

    container "docker.io/thecrick/pipetech_lofreq:amd64-2.1.5"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    lofreq call-parallel --pp-threads $task.cpus -f $fasta -o ${prefix}.vcf $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$( lofreq version 2>&1 | sed -n 's/^version: //p' )
    END_VERSIONS
    """
}

// lofreq call-parallel --pp-threads $task.cpus -f $fasta -o ${prefix}.vcf $bam
// lofreq call -f $fasta -o ${prefix}.vcf $bam