//
// Generate consesnus and call variants from nanopore reads
//

include { MEDAKA_INFERENCE } from '../../../modules/local/medaka/inference/main'
include { MEDAKA_VCF       } from '../../../modules/local/medaka/vcf/main'
include { ARTIC_VCF_MERGE  } from '../../../modules/local/artic/vcf_merge/main'


workflow NANOPORE_VARCALL {
    take:
    trimmed_bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
    primer_trimmed_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    primer_bed             // file
    val_pool_reads         // val
    reference              // channel: [ val(meta), path(fasta), path(fai) ]


    main:
    ch_versions = Channel.empty()

    //
    // CHANNEL: Get primer pool ids and combine with bam
    //
    ch_pool_trimmed_bam = Channel.empty()
    if(val_pool_reads) {
        ch_pool_trimmed_bam = primer_bed
        .splitCsv(sep: '\t')
        .map{[it[0][4]]}
        .unique()
        .flatten()
        .combine(trimmed_bam_bai)
        .map{
            def newMeta = it[1].clone()
            newMeta.pool = it[0]
            [it[0], newMeta, it[2], it[3]]
        }
    } else {
        ch_pool_trimmed_bam = trimmed_bam_bai.map{[[], it[0], it[1]]}
    }

    //
    // MODULE: Generate inference
    //
    MEDAKA_INFERENCE (
        ch_pool_trimmed_bam.map{[it[1], it[2], it[3]]},
        ch_pool_trimmed_bam.map{[it[0]]}
    )
    ch_versions = ch_versions.mix(MEDAKA_INFERENCE.out.versions)
    ch_hdf      = MEDAKA_INFERENCE.out.hdf

    //
    // MODULE: Generate vcf
    //
    MEDAKA_VCF (
        ch_hdf,
        reference.collect()
    )
    ch_versions = ch_versions.mix(MEDAKA_VCF.out.versions)
    ch_vcf      = MEDAKA_VCF.out.vcf

    //
    // MODULE: Merge vcfs
    //
    ch_vcf_pools = ch_vcf
        .map{[it[0].id, it[0], it[1], it[0].pool]}
        .groupTuple()
        .map{
            it[1][0].remove('pool')
            [it[1][0], it[2], it[3]]
        }
    ARTIC_VCF_MERGE (
        ch_vcf_pools,
        primer_bed
    )

    emit:

    versions = ch_versions.ifEmpty(null)
}
