//
// Generate consesnus and call variants from illumina reads
//

include { LOFREQ_CALL                          } from '../../../modules/local/lofreq/call/main'
include { TABIX_BGZIPTABIX as INDEX_LOFREQ_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_PASS_FAIL_SPLIT             } from '../../../modules/local/bcftools/pass_fail_split/main'
include { ARTIC_MAKE_DEPTH_MASK                } from '../../../modules/local/artic/make_depth_mask/main'
include { ARTIC_MASK                           } from '../../../modules/local/artic/mask/main'
include { BCFTOOLS_CONSENSUS                   } from '../../../modules/nf-core/bcftools/consensus/main'
include { LINUX_COMMAND as RENAME_FASTA        } from '../../../modules/local/linux/command/main'
include { FREEBAYES_CALL                       } from '../../../modules/local/freebayes/call/main'
include { TABIX_BGZIPTABIX as INDEX_FREEBAYES  } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow ILLUMINA_VARCALL {
    take:
    bam_bai_fasta_fai // channel: [ val(meta), path(bam), path(bai), path(fasta), path(fai) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Link fasta to bam/bai
    //

    //
    // MODULE: Call low frequency variants
    //
    LOFREQ_CALL (
        bam_bai_fasta_fai.map{[it[0], it[1], it[2]]},
        bam_bai_fasta_fai.map{[it[0], it[3], it[4]]},
    )
    ch_versions   = ch_versions.mix(LOFREQ_CALL.out.versions)
    ch_lofreq_vcf = LOFREQ_CALL.out.vcf

    //
    // MODULE: Gzip and index the lofreq VCF
    //
    INDEX_LOFREQ_VCF (
        ch_lofreq_vcf
    )
    ch_versions       = ch_versions.mix(INDEX_LOFREQ_VCF.out.versions)
    ch_lofreq_vcf_tbi = INDEX_LOFREQ_VCF.out.gz_tbi

    //
    // MODULE: Split VCF into pass/fail files passed on filter
    //
    BCFTOOLS_PASS_FAIL_SPLIT (
        ch_lofreq_vcf
    )
    ch_versions   = ch_versions.mix(BCFTOOLS_PASS_FAIL_SPLIT.out.versions)
    ch_lofreq_vcf = BCFTOOLS_PASS_FAIL_SPLIT.out.pass_vcf

    //
    // MODULE: Make depth mask
    //
    ARTIC_MAKE_DEPTH_MASK (
        bam_bai_fasta_fai.map{[it[0], it[1], it[2]]},
        bam_bai_fasta_fai.map{[it[0], it[3], it[4]]},
    )
    ch_depth_mask = ARTIC_MAKE_DEPTH_MASK.out.mask

    //
    // MODULE: Build a pre-conensus mask of the coverage mask and N's where the failed variants are
    //
    ch_artic_mask = BCFTOOLS_PASS_FAIL_SPLIT.out.fail_vcf
        .map { [it[0].id, it ]}
        .join ( ch_depth_mask.map { [it[0].id, it[1]] })
        .join ( bam_bai_fasta_fai.map { [it[0].id, it[3], it[4]] })
        .map{ [it[1][0], it[1][1], it[2], it[3], it[4]] }
    ARTIC_MASK (
        ch_artic_mask.map{[it[0], it[1]]},
        ch_artic_mask.map{[it[0], it[2]]},
        ch_artic_mask.map{[it[0], it[3], it[4]]},
    )
    ch_preconsensus_mask = ARTIC_MASK.out.fasta

    //
    // MODULE: Call the consensus sequence
    //
    ch_vcf_tbi_fasta_mask = ch_lofreq_vcf_tbi
        .map { [it[0].id, it ]}
        .join ( ch_preconsensus_mask.map { [it[0].id, it[1]] })
        .join ( ch_depth_mask.map { [it[0].id, it[1]] })
        .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3]] }
    BCFTOOLS_CONSENSUS (
        ch_vcf_tbi_fasta_mask
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)
    ch_consensus = BCFTOOLS_CONSENSUS.out.fasta

    //
    // MODULE: Rename consensus fasta
    //
    RENAME_FASTA (
        ch_consensus,
        [],
        true,
        "consensus"
    )
    ch_consensus = RENAME_FASTA.out.file

    //
    // MODULE: Call freebayes variants
    //
    FREEBAYES_CALL (
        bam_bai_fasta_fai.map{[it[0], it[1]]},
        bam_bai_fasta_fai.map{[it[0], it[3], it[4]]},
    )
    ch_versions = ch_versions.mix(FREEBAYES_CALL.out.versions)
    ch_freebayes = FREEBAYES_CALL.out.vcf

    //
    // MODULE: Gzip and index the VCF
    //
    INDEX_FREEBAYES (
        ch_freebayes
    )
    ch_versions       = ch_versions.mix(INDEX_FREEBAYES.out.versions)
    ch_freebayes_vcf_tbi = INDEX_FREEBAYES.out.gz_tbi

    //
    // CHANNEL: Generate merged vcf report channels
    //
    ch_vcf_files = ch_lofreq_vcf.map{[it[0], it[1], "lofreq", 1]}
        .mix(ch_freebayes.map{[it[0], it[1], "freebayes", 2]})
    
    emit:
    versions          = ch_versions.ifEmpty(null)
    consensus         = ch_consensus
    lofreq_vcf_tbi    = ch_lofreq_vcf_tbi
    freebayes_vcf_tbi = ch_freebayes_vcf_tbi
    vcf_files         = ch_vcf_files
}
