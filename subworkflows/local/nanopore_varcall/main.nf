//
// Generate consesnus and call variants from nanopore reads
//




workflow NANOPORE_VARCALL {
    take:
    trimmed_bam           // channel: [ val(meta), path(bam) ]
    primer_trimmed_bam    // channel: [ val(meta), path(bam) ]
    primer_bed            // path
    val_pool_reads


    main:
    ch_versions = Channel.empty()

    //
    // CHANNEL: Get primer pool
    //
    ch_primer_pool = primer_bed
    .splitCsv(sep: '\t')
    .map{[it[0][4]]}
    .unique()

    ch_primer_pool | view

    emit:

    versions = ch_versions.ifEmpty(null)
}
