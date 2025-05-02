process SNAPGENE_TO_GFF {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:0.3.42"

    input:
    tuple val(meta), path(snapgene), val(contig)

    output:
    tuple val(meta), path("*.gff"), emit: gff

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/local/bin/python

    from crick_genome_tools.io.snapgene import convert_to_gff

    convert_to_gff("${snapgene}", "${prefix}.gff", "${contig}", "manual")
    """
}
