process GFF_FLU {
    tag "$meta.id"
    label "process_single"

    container "docker.io/thecrick/pipetech_gfflu:0.0.2"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("gfflu-outdir/*.tsv"), emit: tsv
    tuple val(meta), path("gfflu-outdir/*.faa"), emit: faa
    tuple val(meta), path("gfflu-outdir/*.gbk"), emit: gbk
    tuple val(meta), path("gfflu-outdir/*.miniprot.gff"), emit: miniprot_gff
    tuple val(meta), path("gfflu-outdir/*.snpeff.gff"), emit: snpeff_gff

    script:
    """
    file_name="$fasta"
    base_name="\${file_name%.*}"
    gff_name="\${base_name}.gff"
    new_name="\${base_name}.snpeff.gff"

    gfflu $fasta
    mv gfflu-outdir/\${gff_name} gfflu-outdir/\${new_name}
    """
}
