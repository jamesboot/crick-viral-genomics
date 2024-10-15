//
// Prepare a list of primers against a given reference genome
//

include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN      } from '../../../modules/nf-core/blast/blastn/main'

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
    // MODULE: Blast contigs against viral segments
    //
    BLAST_BLASTN (
        [[id:"primers"], primers_fasta],
        ch_blastdb
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)
    ch_blast    = BLAST_BLASTN.out.txt

    emit:

    versions    = ch_versions                     // channel: [ versions.yml ]
}
