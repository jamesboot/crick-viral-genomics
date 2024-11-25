process LINUX_COMMAND {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(input)
    path input2
    val copy_input
    val output_suffix

    output:
    tuple val(meta), path("*.${output_suffix}.*"), emit: file

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def ext      = task.ext.ext ?: 'txt'
    def post_cmd = task.ext.post_cmd ?: 'echo "NO-ARGS"'
    def pre_cmd  = task.ext.pre_cmd ?: ''
    def copy_cmd = ''
    def target_files = input
    if(copy_input) {
        copy_cmd = "for file in $input; do cp \"\$file\" \"\$file.copy\"; done"
        target_files = input.collect{ file -> "${file}.copy"}.join(' ')
    }
    def main_cmd = "cat $target_files | $post_cmd > ${prefix}.${output_suffix}.${ext}"
    if(task.ext.nofile) {
        main_cmd = "cat $target_files | $post_cmd"
    }

    """
    $copy_cmd
    $pre_cmd
    $main_cmd
    """
}
