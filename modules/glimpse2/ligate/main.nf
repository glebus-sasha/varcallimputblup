process GLIMPSE2_LIGATE {
    tag "$sid"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'

    input:
    tuple val(sid), path(phased_variants)
    tuple val(sid), path(list)

    output:
    tuple val(meta), path("*.bcf"), emit: merged_variants

    script:
    def prefix = "${sid}"
    def suffix = "bcf"
    """

    GLIMPSE2_ligate \
        --input $list \
        --thread $task.cpus \
        --output ${prefix}.${suffix}
    """

    stub:
    def prefix = "${sid}"
    def suffix = "bcf"
    """
    touch ${prefix}.${suffix}
    """
}
