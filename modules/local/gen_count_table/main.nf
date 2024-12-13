process GEN_COUNT_TABLE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:0.3.3"

    input:
    tuple val(meta), path(pileup)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/local/bin/python

    from crick_genome_tools.samtools.utils import count_table_from_pileup

    count_table_from_pileup("${pileup}", "./${prefix}.csv")
    """
}
