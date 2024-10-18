//
// Generate consesnus and call variants from nanopore reads
//

include { MEDAKA_INFERENCE } from '../../../modules/local/medaka/inference/main'
include { MEDAKA_VCF       } from '../../../modules/local/medaka/vcf/main'
include { ARTIC_VCF_MERGE  } from '../../../modules/local/artic/vcf_merge/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { LONGSHOT         } from '../../../modules/local/longshot/main'

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
    ch_merged_vcf = ARTIC_VCF_MERGE.out.vcf

    //
    // MODULE: Gzip and index VCF
    //
    TABIX_BGZIPTABIX (
        ch_merged_vcf
    )
    ch_versions   = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)
    ch_vcf_gz_tbi = TABIX_BGZIPTABIX.out.gz_tbi

    //
    // MODULE: Join VCF with reads
    //
    ch_primer_trimmed_bam_bai_vcf_tbi = primer_trimmed_bam_bai
    .map { [it[0].id, it ]}
    .join ( ch_vcf_gz_tbi.map { [it[0].id, it[1], it[2]] })
    .map{ [it[1][0], it[1][1], it[1][2], it[2], it[3]] }

    //
    // MODULE: Call variants with longshot and provide extra metadata
    //
    LONGSHOT (
        ch_primer_trimmed_bam_bai_vcf_tbi.map{[it[0], it[1], it[2]]},
        ch_primer_trimmed_bam_bai_vcf_tbi.map{[it[0], it[3], it[4]]},
        reference.collect()
    )
    ch_versions = ch_versions.mix(LONGSHOT.out.versions)
    ch_vcf      = LONGSHOT.out.vcf


    emit:

    versions = ch_versions.ifEmpty(null)
}


    // # 8) check and filter the VCFs
    // ## if using strict, run the vcf checker to remove vars present only once in overlap regions (this replaces the original merged vcf from the previous step)
    // if args.strict:
    //     cmds.append("bgzip -f %s.merged.vcf" % (args.sample))
    //     cmds.append("tabix -p vcf %s.merged.vcf.gz" % (args.sample))
    //     cmds.append("artic-tools check_vcf --dropPrimerVars --dropOverlapFails --vcfOut %s.merged.filtered.vcf %s.merged.vcf.gz %s 2> %s.vcfreport.txt" % (args.sample, args.sample, bed, args.sample))
    //     cmds.append("mv %s.merged.filtered.vcf %s.merged.vcf" % (args.sample, args.sample))

    // ## if doing the medaka workflow and longshot required, do it on the merged VCF
    // if args.medaka and not args.no_longshot:
    //     cmds.append("bgzip -f %s.merged.vcf" % (args.sample))
    //     cmds.append("tabix -f -p vcf %s.merged.vcf.gz" % (args.sample))
    //     cmds.append("longshot -P 0 -F -A --no_haps --bam %s.primertrimmed.rg.sorted.bam --ref %s --out %s.merged.vcf --potential_variants %s.merged.vcf.gz" % (args.sample, ref, args.sample, args.sample))
