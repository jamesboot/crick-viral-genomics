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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

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
// check_param_list = [
//     run_dir: params.run_dir,
//     sample_sheet: params.samplesheet
// ]
// for (param in check_param_list) {
//     if (!param.value) {
//         exit 1, "Required parameter not specified: ${param.key}"
//     }
//     else {
//         file(param.value, checkIfExists: true)
//     }
// }

// Check non-manditory input parameters to see if the files exist if they have been specified
// check_param_list = [
//     params.bam
// ]
// for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    // INIT:
    // 
    ch_versions    = Channel.empty()
    ch_samplesheet = Channel.empty()


        // 
        // MODULE: MULTIQC
        // 
        // workflow_summary = multiqc_summary(workflow, params)
        // ch_workflow_summary = Channel.value(workflow_summary)

        // ch_multiqc_files = Channel.empty()
        // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_unique_yml.collect())

        // ch_multiqc_files_all = Channel.empty()
        // ch_multiqc_files_all = ch_multiqc_files_all.mix(ch_multiqc_files)
        // ch_multiqc_files_all = ch_multiqc_files_all.mix(ch_fastqc_zip.collect{it[1]}.ifEmpty([]))
        // ch_multiqc_files_all = ch_multiqc_files_all.mix(PYCOQC_ALL.out.json.collect{it[1]}.ifEmpty([]))
        // ch_multiqc_files_all = ch_multiqc_files_all.mix(NANOPLOT_ALL.out.txt.collect{it[1]}.ifEmpty([]))
        // ch_multiqc_files_all = ch_multiqc_files_all.collect().map{ [[id:"all"], it] }

        // ch_multiqc_files_grouped = ch_grouped_fastqc
        //     .map { [ it[0].id, it[0], it ] }
        //     .join ( PYCOQC_GROUPED.out.json.map{ [ it[0].id, it[1] ] } )
        //     .join ( NANOPLOT_GROUPED.out.txt.map{ [ it[0].id, it[1] ] } )
        //     .map { [ it[1], [ it[2][1].flatten(), it[3], it[4] ].flatten() ] }
        //     .combine(ch_multiqc_files.collect())
        //     .map { [ it[0], [it[1], it[2], it[3], it[4]].flatten() ] }

        // MULTIQC_ALL (
        //     ch_multiqc_files_all,
        //     ch_multiqc_config,
        //     [],
        //     []
        // )
        // multiqc_report_all = MULTIQC_ALL.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
//     dump_parameters(workflow, params)

//     if (params.hook_url) {
//         im_notification(workflow, params, projectDir, runid, summary_params, log)
//     }

//     // if (params.email || params.email_on_fail) {
//     //     NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
//     // }
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
