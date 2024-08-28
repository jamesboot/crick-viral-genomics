process PYTHON_BUILD_REFERENCE_FASTA {
    label 'process_single'

    container "docker.io/thecrick/pipetech_seqsim:dev"

    input:
    tuple val(meta), path(refs)
    tuple val(meta2), path(blast)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    output       = task.ext.output ?: 'top_matches.fasta'
    template 'build_reference_fasta.py'
}
