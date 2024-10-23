process CLAIR3_RUN {
    tag "$meta.id"
    label 'process_high'

    container "docker.io/hkubal/clair3:v1.0.10"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)
    val model
    val platform

    output:
    tuple val(meta), path("*pileup.vcf.gz"), path("*pileup.vcf.gz.tbi")                , emit: pileup_gz_tbi
    tuple val(meta), path("*full_alignment.vcf.gz"), path("*full_alignment.vcf.gz.tbi"), emit: full_alignment_gz_tbi
    tuple val(meta), path("*merge_output.vcf.gz"), path("*merge_output.vcf.gz.tbi")    , emit: merge_output_gz_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /opt/bin/run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${fasta} \\
        --threads=${task.cpus} \\
        --model_path="/opt/models/${model}" \\
        --platform=${platform} \\
        --sample_name=${meta.id} \\
        $args \\
        --output=.

    VCF_FILE="merge_output.vcf.gz"
    TBI_FILE="merge_output.vcf.gz.tbi"
    if [ ! -f "\$TBI_FILE" ]; then
        echo "Index file doesn't exist. Creating it..."
        tabix -p vcf "\$VCF_FILE"
    fi

    mv pileup.vcf.gz ${prefix}.pileup.vcf.gz
    mv full_alignment.vcf.gz ${prefix}.full_alignment.vcf.gz
    mv merge_output.vcf.gz ${prefix}.merge_output.vcf.gz
    mv pileup.vcf.gz.tbi ${prefix}.pileup.vcf.gz.tbi
    mv full_alignment.vcf.gz.tbi ${prefix}.full_alignment.vcf.gz.tbi
    mv merge_output.vcf.gz.tbi ${prefix}.merge_output.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(/opt/bin/run_clair3.sh --version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
