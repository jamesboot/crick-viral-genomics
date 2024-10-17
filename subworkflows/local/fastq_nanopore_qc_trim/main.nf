//
// Read QC and pre-processing for nanopore reads
//

include { TOULLIGQC  } from '../../../modules/local/toulligqc/main'

workflow FASTQ_NANOPORE_QC_TRIM {
    take:
    ch_fastq // channel: [ val(meta), path(reads)  ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run toulligqc on reads
    //
    TOULLIGQC (
        ch_fastq,
        [],
        []
    )
    ch_versions = ch_versions.mix(TOULLIGQC.out.versions)

    emit:

    versions = ch_versions.ifEmpty(null) // channel: [ path(versions.yml) ]
}
