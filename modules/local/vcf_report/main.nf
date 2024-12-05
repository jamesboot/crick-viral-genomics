process VCF_REPORT {
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:0.2.22"

    input:
    tuple val(meta), path(vcfs), val(tools)

    output:
    tuple val(meta), path('*.csv'), emit: csv
    path("*.sh")                  , emit: sh

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = task.ext.output ?: "${prefix}_variant_table.csv"
    """
    #!/usr/local/bin/python

    import os
    from crick_genome_tools.io.vcf import generate_merged_vcf_report

    generate_merged_vcf_report([${vcfs.collect{"\"${it}\""}.join(",")}], ${tools.collect{"\"${it}\""}}, "${output}")

    run_script = '#!/bin/bash\\n\\nSCRIPT_DIR=\$(dirname "\$(realpath "\$0")")\\ndocker run -p 8501:8501 -v "\$SCRIPT_DIR:/data" thecrick/pipetech_vg_report:latest streamlit run crick_genome_tools/reporting/interactive_vcf_report.py -- --data_path /data/${output}'
    rep_path = "${prefix}_run_interative_report.sh"
    with open(rep_path, "w") as f:
        f.write(run_script)
    os.chmod(rep_path, 0o755)
    """
}
