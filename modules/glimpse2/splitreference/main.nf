process GLIMPSE2_SPLITREFERENCE {
    tag "${ref_panel.baseName}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_SPLITREFERENCE"
    errorStrategy 'ignore'   

    input:
    tuple val(chr), path(ref_panel), path(ref_panel_index), val(input_region), val(output_region)

    output:
    tuple val(chr), path("*.bin"), emit: bin_ref

    script:
    def prefix = "${chr}_${output_region.replace(":","_")}"
    """
    GLIMPSE2_split_reference \
        --reference ${ref_panel} \
        --input-region $input_region \
        --output-region $output_region \
        --thread $task.cpus \
        --output ${prefix}
    """

    stub:
    """
    touch ${chr}_${output_region.replace(":","_")}.bin
    """
}