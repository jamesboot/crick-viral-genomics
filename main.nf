#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { params_summary_map } from './modules/local/util/logging/main'
include { summary_log        } from './modules/local/util/logging/main'
include { multiqc_summary    } from './modules/local/util/logging/main'
// include { dump_parameters    } from './modules/local/util/logging/main'
// include { im_notification    } from './modules/local/util/logging/main'
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

params.fasta = get_genome_attribute(params, 'fasta')
params.bwa   = get_genome_attribute(params, 'bwa'  )

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
    fasta: params.fasta
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

include { SEQ_SIMULATOR          } from './modules/local/seq_simulator/main'
include { GUNZIP as GUNZIP_FASTA } from './modules/nf-core/gunzip/main'
include { BWA_INDEX              } from './modules/nf-core/bwa/index/main'
include { SPADES_ASSEMBLE        } from './modules/local/spades/assemble/main'
include { FASTQC                 } from './modules/nf-core/fastqc/main'
include { MULTIQC                } from './modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // INIT
    // 
    ch_versions    = Channel.empty()
    ch_samplesheet = Channel.empty()
    ch_fasta       = file(params.fasta, checkIfExists: true)
    ch_bwa         = Channel.empty()
    if(params.bwa) {
        ch_bwa = file(params.bwa, checkIfExists: true)
    }

    SEQ_SIMULATOR (
        ch_seq_sim_refs.map{[ [id: "${params.seq_sim_profile}_test"], it ]},
        ch_seq_sim_config,
        params.seq_sim_profile,
        params.seq_sim_num_reads
    )
    ch_fastq = SEQ_SIMULATOR.out.fastq

    ch_bwa_index = Channel.empty()
    if(params.run_genome) {
        //
        // MODULE: Uncompress genome fasta file if required
        //
        if (ch_fasta.toString().endsWith(".gz")) {
            ch_fasta    = GUNZIP_FASTA ( [ [id:ch_fasta.baseName], ch_fasta ] ).gunzip
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        }
        else {
            ch_fasta = Channel.from( [ [ [id:ch_fasta.baseName], ch_fasta ] ] )
        }

        //
        // MODULES: Uncompress BWA index or generate if required
        //
        if (params.bwa) {
            if (ch_bwa_index.toString().endsWith(".tar.gz")) {
                ch_bwa_index = UNTAR_BWA ( [[:], ch_bwa_index] ).untar
                ch_versions  = ch_versions.mix( UNTAR_BWA.out.versions )
            } else {
                ch_bwa_index = Channel.of([[:], ch_bwa_index])
            }
        }
        else {
            ch_bwa_index = BWA_INDEX ( ch_fasta ).index
            ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }

    // SPADES_ASSEMBLE (
    //     ch_fastq
    // )

    // FASTQC (
    //     ch_fastq
    // )
    // ch_fastqc_zip = FASTQC.out.zip

    // 
    // MODULE: MULTIQC
    // 
    // workflow_summary = multiqc_summary(workflow, params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect())

    // ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_zip.collect{it[1]}.ifEmpty([]))

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
