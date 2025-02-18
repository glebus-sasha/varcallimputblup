process GLIMPSE2_CHUNK {
    tag "${ref_panel_vcf.baseName}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_CHUNK"

    input:
    tuple val(chr), path(ref_panel_vcf), path(ref_panel_index), path(ref_panel_tsv), path(ref_panel_tsv_index)

    output:
    tuple val(chr), path("${chr}.txt"), emit: chunk_chr

    script:
    """
    GLIMPSE2_chunk \
        --input $ref_panel_vcf \
        --region $chr \
        --sequential \
        --threads 60 \
        --output ${chr}.txt
    """

    stub:
    """
    touch ${chr}.txt
    """
}
