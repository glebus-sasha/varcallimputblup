process GLIMPSE2_LIGATE {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'

    input:
    tuple val(meta), path(input_list), path(input_index)

    output:
    tuple val(meta), path("*vcf.gz"), emit: merged_variants

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    def suffix = "vcf.gz"
    """
    printf "%s\\n" $input_list | tr -d '[],' | sort -V > all_files.txt

    GLIMPSE2_ligate \
        --input all_files.txt \
        --thread $task.cpus \
        --output ${prefix}.${suffix}
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    """
    touch ${prefix}.${suffix}
    """
}
