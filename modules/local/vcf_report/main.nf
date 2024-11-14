process VCF_REPORT {
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:latest"

    input:
    tuple val(meta), path(vcfs), val(tools)

    output:
    tuple val(meta), path('*.csv'), emit: csv

    script:
    output = task.ext.output ?: 'variant_table.csv'
    """
    #!/usr/local/bin/python

    from crick_genome_tools.io.vcf import generate_merged_vcf_report

    generate_merged_vcf_report([${vcfs.collect{"\"${it}\""}.join(",")}], ${tools.collect{"\"${it}\""}}, "${output}")
    """
}
