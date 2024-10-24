process LONGSHOT {
    label 'process_single'
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/longshot:1.0.0--87adca2fab9ee777':
        'community.wave.seqera.io/library/longshot:1.0.0--bd7e35fe51062951' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta3), path(fasta), path(fai)
    tuple val(meta2), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    longshot \\
        $args \\
        --bam $bam \\
        --ref $fasta \\
        --potential_variants $vcf \\
        --out ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longshot: \$( longshot --version 2>&1 | sed 's/Longshot: variant caller (SNVs) for long-read sequencing data //g' )
    END_VERSIONS
    """
}
