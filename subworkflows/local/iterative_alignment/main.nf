include { BWA_INDEX               } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM                 } from '../../../modules/local/bwa/mem/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_CONSENSUS      } from '../../../modules/nf-core/samtools/consensus/main'
include { PYTHON_MASK_FASTA       } from '../../../modules/local/python/mask_fasta/main'

workflow ITERATIVE_ALIGMENT {
    take:
    ref       // [val(meta), path(fasta)]
    reads     // [val(meta), path(fastq)]
    // bwa_args  // val
    // bam       // [val(meta), path(bam)]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Index the latest reference
    //
    BWA_INDEX (
        ref
    )
    // ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
    // ch_bwa_index = BWA_INDEX.out.index

    //
    // MODULE: Align to the reference genome (could be a merged set of segments)
    //
    // BWA_MEM (
    //     reads,
    //     ch_bwa_index,
    //     "-T10 -k 19 -B 4 -O 6"
    // )
    // ch_versions = ch_versions.mix(BWA_MEM.out.versions)
    // ch_bam      = BWA_MEM.out.bam

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    // BAM_SORT_STATS_SAMTOOLS (
    //     ch_bam,
    //     ref
    // )
    // ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)
    // ch_bam      = BAM_SORT_STATS_SAMTOOLS.out.bam
    // ch_bai      = BAM_SORT_STATS_SAMTOOLS.out.bai

    //
    // MODULE: Build a consensus
    //
    // SAMTOOLS_CONSENSUS (
    //     ch_bam
    // )
    // ch_versions  = ch_versions.mix(SAMTOOLS_CONSENSUS.out.versions)
    // ch_consensus = SAMTOOLS_CONSENSUS.out.fasta

    //
    // MODULE: Mask against reference
    //
    // PYTHON_MASK_FASTA (
    //     ch_consensus,
    //     ref
    // )
    // ch_versions = ch_versions.mix(PYTHON_MASK_FASTA.out.versions)


    emit:
    ref
    reads
    // bwa_args
    // bam = ch_bam.collect()
    // bai = ch_bai.collect()
    // stats    = BAM_SORT_STATS_SAMTOOLS.out.stats.collect()
    // flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat.collect()
    // idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect()
}
