//
// Generate consesnus and call variants from nanopore reads
//

include { MEDAKA_INFERENCE                     } from '../../../modules/local/medaka/inference/main'
include { MEDAKA_VCF                           } from '../../../modules/local/medaka/vcf/main'
include { MEDAKA_ANNOTATE                      } from '../../../modules/local/medaka/annotate/main'
include { ARTIC_VCF_MERGE                      } from '../../../modules/local/artic/vcf_merge/main'
include { BCFTOOLS_PASS_FAIL_SPLIT             } from '../../../modules/local/bcftools/pass_fail_split/main'
include { TABIX_BGZIPTABIX as INDEX_CONSEN_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { ARTIC_MAKE_DEPTH_MASK                } from '../../../modules/local/artic/make_depth_mask/main'
include { ARTIC_MASK                           } from '../../../modules/local/artic/mask/main'
include { BCFTOOLS_CONSENSUS                   } from '../../../modules/nf-core/bcftools/consensus/main'
include { CLAIR3_RUN                           } from '../../../modules/local/clair3/main'
include { GUNZIP as GUNZIP_VCF                 } from '../../../modules/nf-core/gunzip/main'
include { LOFREQ_CALL                          } from '../../../modules/local/lofreq/call/main'

workflow NANOPORE_VARCALL {
    take:
    trimmed_bam_bai        // channel: [ val(meta), path(bam), path(bai) ]
    primer_trimmed_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    primer_bed             // file
    val_pool_reads         // val
    reference              // channel: [ val(meta), path(fasta), path(fai) ]
    gff                    // file
    clair3_model           // val
    clair3_platform        // val

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
    ch_versions   = ch_versions.mix(MEDAKA_INFERENCE.out.versions)
    ch_medaka_hdf = MEDAKA_INFERENCE.out.hdf

    //
    // MODULE: Generate vcf
    //
    MEDAKA_VCF (
        ch_medaka_hdf,
        reference.collect()
    )
    ch_versions   = ch_versions.mix(MEDAKA_VCF.out.versions)
    ch_medaka_vcf = MEDAKA_VCF.out.vcf

    //
    // MODULE: Join pooled VCF with reads
    //
    ch_pool_trimmed_bam_vcf = ch_pool_trimmed_bam
    .map { [it[1].id + "_" + it[1].pool, it ]}
    .join ( ch_medaka_vcf.map { [it[0].id + "_" + it[0].pool, it[1]] })
    .map{ [it[1][0], it[1][1], it[1][2], it[1][3], it[2]] }

    //
    // MODULE: Generate vcf
    //
    MEDAKA_ANNOTATE (
        ch_pool_trimmed_bam_vcf.map{[it[1], it[4]]},
        ch_pool_trimmed_bam_vcf.map{[it[1], it[2], it[3]]},
        reference.collect(),
        ch_pool_trimmed_bam_vcf.map{[it[0]]}
    )
    ch_versions   = ch_versions.mix(MEDAKA_ANNOTATE.out.versions)
    ch_medaka_vcf = MEDAKA_ANNOTATE.out.vcf

    //
    // MODULE: Merge medaka vcfs
    //
    ch_vcf_pools = ch_medaka_vcf
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
    ch_medaka_vcf = ARTIC_VCF_MERGE.out.vcf

    //
    // MODULE: Split VCF into pass/fail files passed on filter
    //
    BCFTOOLS_PASS_FAIL_SPLIT (
        ch_medaka_vcf
    )
    ch_versions   = ch_versions.mix(BCFTOOLS_PASS_FAIL_SPLIT.out.versions)
    ch_medaka_vcf = BCFTOOLS_PASS_FAIL_SPLIT.out.pass_vcf

    //
    // MODULE: Gzip and index the consensus VCF
    //
    INDEX_CONSEN_VCF (
        ch_medaka_vcf
    )
    ch_versions          = ch_versions.mix(INDEX_CONSEN_VCF.out.versions)
    ch_medaka_vcf_gz_tbi = INDEX_CONSEN_VCF.out.gz_tbi

    //
    // MODULE: Make depth mask
    //
    ARTIC_MAKE_DEPTH_MASK (
        primer_trimmed_bam_bai,
        reference.collect()
    )
    ch_depth_mask = ARTIC_MAKE_DEPTH_MASK.out.mask

    //
    // MODULE: Build a pre-conensus mask of the coverage mask and N's where the failed variants are
    //
    ARTIC_MASK (
        BCFTOOLS_PASS_FAIL_SPLIT.out.fail_vcf,
        ch_depth_mask,
        reference.collect()
    )
    ch_preconsensus_mask = ARTIC_MASK.out.fasta

    //
    // MODULE: Call the consensus sequence
    //
    ch_vcf_tbi_fasta_mask = ch_medaka_vcf_gz_tbi
        .map { [it[0].id, it ]}
        .join ( ch_preconsensus_mask.map { [it[0].id, it[1]] })
        .join ( ch_depth_mask.map { [it[0].id, it[1]] })
        .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3]] }
    BCFTOOLS_CONSENSUS (
        ch_vcf_tbi_fasta_mask
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)
    ch_consensus = BCFTOOLS_CONSENSUS.out.fasta

    // TODO Adjust header with proper fasta header

    //
    // MODULE: Run clair3 variant caller for more accurate variant calling
    //
    CLAIR3_RUN (
        primer_trimmed_bam_bai,
        reference.collect(),
        clair3_model,
        clair3_platform
    )
    ch_versions          = ch_versions.mix(CLAIR3_RUN.out.versions)
    ch_clair3_vcf_gz_tbi = CLAIR3_RUN.out.merge_output_gz_tbi
    ch_clair3_vcf_unzip  = CLAIR3_RUN.out.pileup_gz_tbi
                            .mix(CLAIR3_RUN.out.full_alignment_gz_tbi)
                            .mix(CLAIR3_RUN.out.merge_output_gz_tbi)

    //
    // MODULE: Unzip clair3 VCF files
    //
    GUNZIP_VCF (
        ch_clair3_vcf_unzip.map{[it[0], it[1]]}
    )
    ch_versions = ch_versions.mix(GUNZIP_VCF.out.versions)

    //
    // MODULE: Call low frequency variants
    //
    LOFREQ_CALL (
        primer_trimmed_bam_bai,
        reference.collect(),
    )
    ch_versions   = ch_versions.mix(LOFREQ_CALL.out.versions)
    ch_lofreq_vcf = LOFREQ_CALL.out.vcf

    emit:
    versions       = ch_versions.ifEmpty(null)
    consensus      = ch_consensus
    clair3_vcf_tbi = ch_clair3_vcf_gz_tbi
    lofreq_vcf     = ch_lofreq_vcf
    medaka_vcf_tbi = ch_medaka_vcf_gz_tbi
}
