process PYTHON_MASK_FASTA {
    label 'process_single'

    container "docker.io/thecrick/pipetech_seqsim:dev"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(ref)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    output       = task.ext.output ?: 'masked.fasta'
    template 'mask_fasta.py'
}
