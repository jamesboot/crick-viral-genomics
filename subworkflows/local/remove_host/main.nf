//
// Remove host from illumina reads
//

include { GUNZIP as GUNZIP_FASTA                         } from '../../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_BWA_INDEX                       } from '../../../modules/nf-core/untar/main'
include { BWA_INDEX as BWA_INDEX_HOST                    } from '../../../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_ALIGN_HOST                      } from '../../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_INDEX                                 } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                                 } from '../../../modules/nf-core/minimap2/align/main'
include { BAM_SORT_STATS_SAMTOOLS as BAM_HOST_SORT_STATS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HOST            } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_VIRAL           } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_VIRAL         } from '../../../modules/nf-core/samtools/fastq/main'

workflow REMOVE_HOST {
    take:
    ch_fastq        // channel: [ val(meta), path(reads) ]
    host_fasta      // file: [ path(fasta) ]
    host_bwa_index  // file: [ path(index) ]
    mode            // val: ont/illumina

    main:
    ch_versions = Channel.empty()

    ch_host_fasta = Channel.empty()
    if(host_fasta) {
        //
        // MODULE: Uncompress host genome fasta file if required
        //
        if (host_fasta.toString().endsWith(".gz")) {
            ch_host_fasta = GUNZIP_FASTA ( [ [id:host_fasta.baseName], host_fasta ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        }
        else {
            ch_host_fasta = Channel.from( [ [ [id:host_fasta.baseName], host_fasta ] ] )
        }
    }

    //
    // MODULES: Uncompress BWA index or generate if required
    //
    ch_host_bwa_index = Channel.empty()
    if (mode == "illumina") {
        if (host_bwa_index) {
            if (host_bwa_index.toString().endsWith(".tar.gz")) {
                UNTAR_BWA_INDEX ( [[:], host_bwa_index] )
                ch_versions       = ch_versions.mix( UNTAR_BWA_INDEX.out.versions )
                ch_host_bwa_index = UNTAR_BWA_INDEX.out.untar

            } else {
                ch_host_bwa_index = Channel.of([[:], host_bwa_index])
            }
        }
        else {
            BWA_INDEX_HOST (ch_host_fasta )
            ch_versions       = ch_versions.mix(BWA_INDEX_HOST.out.versions)
            ch_host_bwa_index = BWA_INDEX_HOST.out.index
        }

        //
        // MODULE: Map raw reads to host
        //
        BWA_ALIGN_HOST (
            ch_fastq,
            ch_host_bwa_index.collect(),
            ch_host_fasta.collect(),
            "view"
        )
        ch_versions = ch_versions.mix(BWA_ALIGN_HOST.out.versions)
        ch_host_bam = BWA_ALIGN_HOST.out.bam
    } else if(mode == "ont") {
        //
        // MODULE: Build minimap index
        //
        MINIMAP2_INDEX ( ch_host_fasta )
        ch_versions      = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_host_mm_index = MINIMAP2_INDEX.out.index

        //
        // MODULE: Map raw reads to host
        //
        MINIMAP2_ALIGN (
            ch_fastq,
            ch_host_mm_index.collect(),
            true,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
        ch_host_bam = MINIMAP2_ALIGN.out.bam
    }

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_HOST_SORT_STATS (
        ch_host_bam,
        ch_host_fasta.collect()
    )
    ch_versions          = ch_versions.mix(BAM_HOST_SORT_STATS.out.versions)
    ch_host_bam          = BAM_HOST_SORT_STATS.out.bam
    ch_host_bai          = BAM_HOST_SORT_STATS.out.bai
    ch_host_bam_stats    = BAM_HOST_SORT_STATS.out.stats
    ch_host_bam_flagstat = BAM_HOST_SORT_STATS.out.flagstat
    ch_host_bam_idxstats = BAM_HOST_SORT_STATS.out.idxstats

    //
    // CHANNEL: Combine bam and bai files
    //
    ch_host_bam_bai = ch_host_bam
    .map { row -> [row[0].id, row ].flatten()}
    .join ( ch_host_bai.map { row -> [row[0].id, row ].flatten()} )
    .map { row -> [row[1], row[2], row[4]] }

    //
    // MODULE: Filter for unmapped reads
    //
    SAMTOOLS_VIEW_HOST (
        ch_host_bam_bai,
        [[],[]],
        []
    )
    ch_versions  = ch_versions.mix(SAMTOOLS_VIEW_HOST.out.versions)
    ch_viral_bam = SAMTOOLS_VIEW_HOST.out.bam

    //
    // MODULE: Sort by name
    //
    SAMTOOLS_SORT_VIRAL (
        ch_viral_bam,
        [[],[]]
    )
    ch_versions  = ch_versions.mix(SAMTOOLS_VIEW_HOST.out.versions)
    ch_viral_bam = SAMTOOLS_SORT_VIRAL.out.bam

    //
    // MODULE: Convert to fastq
    //
    SAMTOOLS_FASTQ_VIRAL (
        ch_viral_bam,
        false
    )
    ch_versions    = ch_versions.mix(SAMTOOLS_FASTQ_VIRAL.out.versions)

    // Reads come out of different channels depending if we have paried end or single-end
    if(mode == "illumina") {
        ch_viral_fastq = SAMTOOLS_FASTQ_VIRAL.out.fastq
    } else {
        ch_viral_fastq = SAMTOOLS_FASTQ_VIRAL.out.other
    }

    emit:
    viral_fastq       = ch_viral_fastq
    host_bam          = ch_host_bam
    host_bai          = ch_host_bai
    host_bam_stats    = ch_host_bam_stats
    host_bam_flagstat = ch_host_bam_flagstat
    host_bam_idxstats = ch_host_bam_idxstats
    versions          = ch_versions
}
