process GLIMPSE2_SPLITREFERENCE {
    tag "${ref_panel.baseName}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_SPLITREFERENCE"

    input:
    path(ref_panel)
    path(ref_panel_index)
    path(chunk_chr)
    val(input_region)
    val(output_region)

    output:
        tuple val(meta), path("*.bin"), emit: bin_ref
        path "versions.yml"           , emit: versions

    script:
    """
    GLIMPSE2_split_reference --reference ${ref_panel} --input-region  chr1:154644662-1248956422 --output-region   chr1:155644650-1248956422 --output reference_panel/split/1000GP.chr22.noNA12878
    """

    stub:
    """
        GLIMPSE2_split_reference \
        --reference $ref_panel \
        --input-region $input_region \
        --output-region $output_region \
        --thread $task.cpus \
        --output ${prefix}
    touch ${ref_panel.baseName}.bin
    """
}