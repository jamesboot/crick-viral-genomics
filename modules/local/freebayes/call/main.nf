process FREEBAYES_CALL {
    tag "$meta.id"
    label "process_single"

    container "community.wave.seqera.io/library/freebayes:1.3.8--5fb5db8f54462c80"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    freebayes \
      -f $fasta \
      $args \
      $bam > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$( freebayes --version 2>&1 | sed -n 's/^version:  \\(v[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\)\$/\\1/p' )
    END_VERSIONS
    """
}
