process GLIMPSE2_PHASE {
    tag ""
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_PHASE"
//    errorStrategy 'ignore'   

    input:
        path ref_panel_bin
        path ref_panel_index
        tuple (file), path(bam), path(bamindex)

    output:
        path("*.{vcf,vcf.gz,bcf,bgen}"), emit: phased_variants

    script:
    """
    GLIMPSE2_phase \
        --reference $ref_panel_bin \
        --bam-file $bam \
        --thread $task.cpus \
        --output "${ref_panel_bin.baseName}.bcf"

    """

    stub:
    """
    touch ${prefix}.${suffix}
    """
}
