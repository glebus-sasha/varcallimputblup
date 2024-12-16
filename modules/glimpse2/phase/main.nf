process GLIMPSE2_PHASE {
    tag "${ref_panel_bin.baseName}_${bam.baseName}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_PHASE"
//    errorStrategy 'ignore'   

    input:
        tuple val(sid), path(bam), path(bamindex), path(ref_panel_bin), val(chr), path(ref_panel_index)

    output:
        tuple val(sid), path("*.bcf"), path("*.bcf.csi"), emit: phased_variants

    script:
    """
    GLIMPSE2_phase \
        --reference $ref_panel_bin \
        --bam-file $bam \
        --thread $task.cpus \
        --output "${sid}_${ref_panel_bin.baseName}.bcf"
    """

    stub:
    """
    touch ${sid}_${ref_panel_bin.baseName}.bcf
    touch "${sid}.txt"
    """
}
