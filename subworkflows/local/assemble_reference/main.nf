//
// Assemble reference genome from viral fasta files
//

include { SPADES_ASSEMBLE       } from '../../../modules/local/spades/assemble/main'
include { CDHIT_CDHIT           } from '../../../modules/nf-core/cdhit/cdhit/main'
include { SEQTK_FASTA           } from '../../../modules/local/seqtk/fasta/main'
include { BLAST_MAKEBLASTDB     } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN          } from '../../../modules/nf-core/blast/blastn/main'
include { BUILD_REFERENCE_FASTA } from '../../../modules/local/build_reference_fasta/main'


workflow ASSEMBLE_REFERENCE {
    take:
    fastq       // channel: [ val(meta), path(reads) ]
    viral_fasta // channel: [ val(meta), path(reads) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: First assembly of non-host reads
    //
    SPADES_ASSEMBLE (
        fastq
    )
    ch_contigs = SPADES_ASSEMBLE.out.contigs

    //
    // MODULE: Cluster contigs
    //
    CDHIT_CDHIT (
        ch_contigs
    )
    ch_versions          = ch_versions.mix(CDHIT_CDHIT.out.versions)
    ch_clustered_contigs = CDHIT_CDHIT.out.fasta

    //
    // MODULE: Filter small contigs
    //
    SEQTK_FASTA (
        ch_clustered_contigs
    )
    ch_versions         = ch_versions.mix(SEQTK_FASTA.out.versions)
    ch_filtered_contigs = SEQTK_FASTA.out.fasta

    //
    // MODULE: Build blastdb of viral segments
    //
    BLAST_MAKEBLASTDB (
        viral_fasta
    )
    ch_versions      = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
    ch_viral_blastdb = BLAST_MAKEBLASTDB.out.db

    //
    // MODULE: Blast contigs against viral segments
    //
    BLAST_BLASTN (
        ch_filtered_contigs,
        ch_viral_blastdb.collect()
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)
    ch_blast    = BLAST_BLASTN.out.txt

    //
    // MODULE: Build reference fasta from top blast hits
    //
    BUILD_REFERENCE_FASTA (
        viral_fasta,
        ch_blast
    )
    ch_viral_ref = BUILD_REFERENCE_FASTA.out.fasta

    //
    // CHANNEL: Join ref to reads
    //
    ch_fastq_fasta = fastq
    .map { [it[0].id, it ]}
    .join ( ch_viral_ref.map { [it[1].simpleName, it[1]] })
    .map{ [it[1][0], it[1][1], it[2]] }

    //
    // CHANNEL: Join ref to metadata
    //
    ch_viral_ref = fastq
    .map { [it[0].id, it[0] ]}
    .join ( ch_viral_ref.map { [it[1].simpleName, it[1]] })
    .map{ [it[1], it[2]] }

    emit:
    viral_ref   = ch_viral_ref
    fastq_fasta = ch_fastq_fasta
    versions = ch_versions.ifEmpty(null)
}
