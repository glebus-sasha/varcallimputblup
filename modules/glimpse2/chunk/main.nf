process GLIMPSE2_CHUNK {
    tag "${ref_panel.baseName}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_CHUNK"

    input:
    tuple val(chr), path(ref_panel), path(ref_panel_index)

    output:
    path("${ref_panel.simpleName}.txt"), emit: chunk_chr

    script:
    """
    GLIMPSE2_chunk \
        --input $ref_panel \
        --region $chr \
        --sequential \
        --threads 60 \
        --output ${ref_panel.simpleName}.txt
    """

    stub:
    """
    touch ${ref_panel.simpleName}.txt
    """
}
