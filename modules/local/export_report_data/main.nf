process EXPORT_REPORT_DATA {
    tag "$run_id"
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:0.3.22"

    input:
    val(run_id)
    path("data/toulligqc/*")
    path("data/samtools_host/*")
    path("data/samtools_contaminent/*")

    output:
    path("*.pkl"), emit: pkl

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/local/bin/python

    from crick_genome_tools.reporting.report_data_parser import ReportDataParser

    parser = ReportDataParser("./data")
    parser.get_data()
    parser.save_data("${run_id}.pkl")
    """
}
