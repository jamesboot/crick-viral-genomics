process BUILD_REFERENCE_FASTA {
    label 'process_single'

    container "docker.io/thecrick/pipetech_genome_tools:latest"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(blast)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    output = task.ext.output ?: 'top_matches.fasta'
    """
    #!/usr/local/bin/python

    from crick_genome_tools.blast.utils import build_fasta_from_top_hits
    from crick_genome_tools.io.fasta import Fasta

    fasta_data = build_fasta_from_top_hits("${fasta}", "${blast}")
    Fasta.write_fasta_file(fasta_data, "${output}")
    """
}
