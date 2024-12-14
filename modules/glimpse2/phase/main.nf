process GLIMPSE2_PHASE {
    tag "${ref_panel_bin.baseName}.${bam.baseName}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_PHASE"
//    errorStrategy 'ignore'   

    input:
        each path(ref_panel_bin)
        tuple val(chunk_chr), path(ref_panel_index)
        tuple val(sid), path(bam), path(bamindex)

    output:
        tuple val(sid), path("*.bcf"), emit: phased_variants

    script:
    def region = "${output_region.replace(":","_")}"
    def prefix = "${meta.id}_${region}"
    def suffix = "bcf"
    """
    GLIMPSE2_phase \
        --reference $ref_panel_bin \
        --bam-file $bam \
        --thread $task.cpus \
        --output "${sid}_${chunk_chr}.bcf"

    """

    stub:
    """
    touch ${prefix}.${suffix}
    """
}
