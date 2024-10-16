//
// Read QC and trimming for nanopore reads
//

// include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
// include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
// include { FASTP                 } from '../../../modules/nf-core/fastp/main'


workflow FASTQ_NANOPORE_QC_TRIM {
    take:
    ch_reads              // channel: [ val(meta), path(reads)  ]
    // ch_adapter_fasta      // channel: [ path(fasta) ]
    // val_save_trimmed_fail // value: boolean
    // val_save_merged       // value: boolean
    // val_skip_fastp        // value: boolean
    // val_skip_fastqc       // value: boolean

    main:

    ch_versions = Channel.empty()


    emit:
    // reads             = ch_trim_reads         // channel: [ val(meta), path(reads) ]
    // trim_json         = ch_trim_json          // channel: [ val(meta), path(json) ]
    // trim_html         = ch_trim_html          // channel: [ val(meta), path(html) ]
    // trim_log          = ch_trim_log           // channel: [ val(meta), path(log) ]
    // trim_reads_fail   = ch_trim_reads_fail    // channel: [ val(meta), path(fastq.gz) ]
    // trim_reads_merged = ch_trim_reads_merged  // channel: [ val(meta), path(fastq.gz) ]

    // fastqc_raw_html  = ch_fastqc_raw_html    // channel: [ val(meta), path(html) ]
    // fastqc_raw_zip   = ch_fastqc_raw_zip     // channel: [ val(meta), path(zip) ]
    // fastqc_trim_html = ch_fastqc_trim_html   // channel: [ val(meta), path(html) ]
    // fastqc_trim_zip  = ch_fastqc_trim_zip    // channel: [ val(meta), path(zip) ]

    versions = ch_versions.ifEmpty(null) // channel: [ path(versions.yml) ]
}
