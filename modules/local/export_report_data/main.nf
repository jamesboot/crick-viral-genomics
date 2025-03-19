process EXPORT_REPORT_DATA {
    tag "$run_id"
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:0.3.22"

    input:
    val(run_id)
    val(summary)
    path("data/samplesheet/*")
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

    import json
    from crick_genome_tools.reporting.report_data_parser import ReportDataParser

    summary_dict = json.loads('${summary}')
    with open("./data/summary.json", "w") as f:
        json.dump(summary_dict, f)

    parser = ReportDataParser("./data")
    parser.get_data()
    parser.save_data("${run_id}.pkl")
    """
}
