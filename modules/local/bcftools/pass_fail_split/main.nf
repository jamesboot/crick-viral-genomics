process BCFTOOLS_PASS_FAIL_SPLIT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.filtered.vcf'), emit: vcf
    tuple val(meta), path('*.pass.vcf')    , emit: pass_vcf
    tuple val(meta), path('*.fail.vcf')    , emit: fail_vcf
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filter = task.ext.filter ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools filter \\
        -e '$filter' \\
        -s FAIL \\
        -S . \\
        $vcf \\
        > ${prefix}.filtered.vcf

    bcftools view -f PASS ${prefix}.filtered.vcf > ${prefix}.pass.vcf
    bcftools view -f FAIL ${prefix}.filtered.vcf > ${prefix}.fail.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
