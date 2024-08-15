// include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
// include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
// include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow ITERATIVE_ALIGMENT {
    take:
   

    main:

    ch_versions = Channel.empty()


    emit:
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    // stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    // flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    // idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
