//
// Prepare a list of primers against a given reference genome
//

include { BLAST_MAKEBLASTDB            } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                 } from '../../../modules/nf-core/blast/blastn/main'
include { LINUX_COMMAND as BLAST_2_BED }  from '../../../modules/local/linux/command/main'

workflow PREPARE_PRIMERS {
    take:
    ref_fasta     // channel: [ val(meta), [ fasta ] ]
    primers_fasta // path   : primer fasta file
    primers_csv   // path   : primer csv file

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Build blastdb of viral segments
    //
    BLAST_MAKEBLASTDB (
        ref_fasta
    )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
    ch_blastdb  = BLAST_MAKEBLASTDB.out.db

    //
    // MODULE: Blast primers against ref
    //
    BLAST_BLASTN (
        [[id:"primers"], primers_fasta],
        ch_blastdb
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)
    ch_blast    = BLAST_BLASTN.out.txt

    //
    // CHANNEL: Add sample meta to blast results
    //
    ch_blast_2_bed = ch_blast
        .map{ [ [id:it[1].toString().split('/')[-1].split('\\.')[0]], it[1] ]}

    //
    // MODULE: Convert to bed
    //
    BLAST_2_BED (
        ch_blast_2_bed,
        [],
        false,
        "primers"
    )
    ch_primer_bed = BLAST_2_BED.out.file

    emit:
    versions   = ch_versions
    primer_bed = ch_primer_bed
}
