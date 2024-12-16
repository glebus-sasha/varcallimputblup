process GLIMPSE2_LIGATE {
    tag "$sid"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_LIGATE"

    input:
    tuple val(sid), path(phased_variants), path(phased_variants_index)
    
    output:
    tuple val(sid), path("*.bcf"), emit: merged_variants

    script:
    def prefix = "${sid}"
    def suffix = "_imputed.bcf"
    """
    ls -1v *.bcf > files.txt
    GLIMPSE2_ligate \
        --input files.txt \
        --thread $task.cpus \
        --output ${sid}${suffix}
    """

    stub:
    def prefix = "${sid}"
    def suffix = "bcf"
    """
    touch ${prefix}.${suffix}
    """
}
