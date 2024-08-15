#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
// nextflow.preview.recursion = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { params_summary_map } from './modules/local/util/logging/main'
include { summary_log        } from './modules/local/util/logging/main'
include { multiqc_summary    } from './modules/local/util/logging/main'
include { get_genome_attribute } from './modules/local/util/references/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_seq_sim_config = file(params.seq_sim_config, checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.host_fasta = get_genome_attribute(params, 'fasta')
params.host_bwa   = get_genome_attribute(params, 'bwa'  )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info summary_log(workflow, params, params.debug, params.monochrome_logs)
def summary_params = params_summary_map(workflow, params, params.debug)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    host_fasta: params.host_fasta
]
for (param in check_param_list) {
    if (!param.value) {
        exit 1, "Required parameter not specified: ${param.key}"
    }
    else {
        file(param.value, checkIfExists: true)
    }
}

// Check non-manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    params.sample_sheet,
    params.seq_sim_ref_dir,
    params.seq_sim_config
]
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }


// Load seq sim ref dir
ch_seq_sim_refs = Channel.empty()
if (params.generate_reads) {
    ch_seq_sim_refs = Channel.from(file(params.seq_sim_ref_dir, checkIfExists: true))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQ_SIMULATOR                        } from './modules/local/seq_simulator/main'
include { GUNZIP as GUNZIP_FASTA               } from './modules/nf-core/gunzip/main'
include { FASTQC                               } from './modules/nf-core/fastqc/main'
include { BWA_INDEX                            } from './modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_ALIGN_HOST            } from './modules/nf-core/bwa/mem/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_HOST  } from './modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_VIRAL } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FASTQ                       } from './modules/nf-core/samtools/fastq/main'
include { SPADES_ASSEMBLE                      } from './modules/local/spades/assemble/main'
include { LINUX_COMMAND as MERGE_REFS          } from './modules/local/linux/command/main'
include { BLAST_MAKEBLASTDB                    } from './modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                         } from './modules/nf-core/blast/blastn/main'
include { PYTHON_BUILD_REFERENCE_FASTA         } from './modules/local/python/build_reference_fasta/main'
include { MULTIQC                              } from './modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_SORT_STATS_SAMTOOLS as BAM_HOST_SORT_STATS } from './subworkflows/nf-core/bam_sort_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // INIT
    //
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_samplesheet   = Channel.empty()
    ch_host_fasta    = file(params.host_fasta, checkIfExists: true)
    ch_host_bwa      = Channel.empty()
    if(params.host_bwa) {
        ch_host_bwa = file(params.host_bwa, checkIfExists: true)
    }
    ch_viral_fasta = Channel.fromPath(params.viral_fasta).toSortedList().map{[[id:"fasta"], it]}

    //
    // MODULE: Uncompress genome fasta file if required
    //
    SEQ_SIMULATOR (
        ch_seq_sim_refs.map{[ [id: "${params.seq_sim_profile}_test"], it ]},
        ch_seq_sim_config,
        params.seq_sim_profile,
        params.seq_sim_num_reads
    )
    ch_fastq = SEQ_SIMULATOR.out.fastq

    ch_host_bwa_index = Channel.empty()
    if(params.run_genome) {
        //
        // MODULE: Uncompress genome fasta file if required
        //
        if (ch_host_fasta.toString().endsWith(".gz")) {
            ch_host_fasta = GUNZIP_FASTA ( [ [id:ch_host_fasta.baseName], ch_host_fasta ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        }
        else {
            ch_host_fasta = Channel.from( [ [ [id:ch_host_fasta.baseName], ch_host_fasta ] ] )
        }

        //
        // MODULES: Uncompress BWA index or generate if required
        //
        if (params.host_bwa) {
            if (ch_host_bwa_index.toString().endsWith(".tar.gz")) {
                UNTAR_BWA ( [[:], ch_host_bwa_index] )
                ch_versions       = ch_versions.mix( UNTAR_BWA.out.versions )
                ch_host_bwa_index = UNTAR_BWA.out.untar

            } else {
                ch_host_bwa_index = Channel.of([[:], ch_host_bwa_index])
            }
        }
        else {
            BWA_INDEX ( ch_host_fasta )
            ch_versions       = ch_versions.mix(BWA_INDEX.out.versions)
            ch_host_bwa_index = BWA_INDEX.out.index
        }
    }

    //
    // MODULE: Fastqc on raw reads
    //
    FASTQC (
        ch_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

    //
    // MODULE: Map raw reads to host
    //
    BWA_ALIGN_HOST (
        ch_fastq,
        ch_host_bwa_index,
        ch_host_fasta,
        "view"
    )
    ch_versions = ch_versions.mix(BWA_ALIGN_HOST.out.versions)
    ch_host_bam = BWA_ALIGN_HOST.out.bam

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_HOST_SORT_STATS (
        ch_host_bam,
        ch_host_fasta
    )
    ch_versions           = ch_versions.mix(BAM_HOST_SORT_STATS.out.versions)
    ch_host_bam           = BAM_HOST_SORT_STATS.out.bam
    ch_host_bai           = BAM_HOST_SORT_STATS.out.bai
    ch_host_bam_stats     = BAM_HOST_SORT_STATS.out.stats
    ch_host_bams_flagstat = BAM_HOST_SORT_STATS.out.flagstat
    ch_host_bam_idxstats  = BAM_HOST_SORT_STATS.out.idxstats
    ch_multiqc_files      = ch_multiqc_files.mix(ch_host_bam_stats.collect{it[1]})
    ch_multiqc_files      = ch_multiqc_files.mix(ch_host_bams_flagstat.collect{it[1]})
    ch_multiqc_files      = ch_multiqc_files.mix(ch_host_bam_idxstats.collect{it[1]})

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
        ch_host_fasta,
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_HOST.out.versions)
    ch_bam = SAMTOOLS_VIEW_HOST.out.bam

    //
    // MODULE: Sort by name
    //
    SAMTOOLS_SORT_VIRAL (
        ch_bam,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_HOST.out.versions)
    ch_bam = SAMTOOLS_SORT_VIRAL.out.bam

    //
    // MODULE: Convert to fastq
    //
    SAMTOOLS_FASTQ (
        ch_bam,
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
    ch_fastq = SAMTOOLS_FASTQ.out.fastq

    //
    // MODULE: First assembly of non-host reads
    //
    SPADES_ASSEMBLE (
        ch_fastq
    )
    // ch_versions = ch_versions.mix(SPADES_ASSEMBLE.out.versions)
    ch_contigs = SPADES_ASSEMBLE.out.contigs

    //
    // MODULE: Concat the reference files into one file
    //
    MERGE_REFS (
        ch_viral_fasta,
        [],
        true
    )
    ch_merged_viral_fasta = MERGE_REFS.out.file

    //
    // MODULE: Build blastdb of viral segments
    //
    BLAST_MAKEBLASTDB (
        ch_merged_viral_fasta
    )
    ch_versions      = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
    ch_viral_blastdb = BLAST_MAKEBLASTDB.out.db

    //
    // MODULE: Blast contigs against viral segments
    //
    BLAST_BLASTN (
        ch_contigs,
        ch_viral_blastdb
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)
    ch_blast    = BLAST_BLASTN.out.txt

    //
    // MODULE: Build reference fasta from top blast hits
    //
    PYTHON_BUILD_REFERENCE_FASTA (
        ch_merged_viral_fasta,
        ch_blast
    )
    ch_versions       = ch_versions.mix(PYTHON_BUILD_REFERENCE_FASTA.out.versions)
    ch_top_hits_fasta = PYTHON_BUILD_REFERENCE_FASTA.out.fasta


    // 
    // MODULE: MULTIQC
    // 
    // workflow_summary = multiqc_summary(workflow, params)
    // ch_workflow_summary = Channel.value(workflow_summary)
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect())

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config,
    //     [],
    //     []
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
