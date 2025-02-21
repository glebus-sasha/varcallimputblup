process GLIMPSE2_CONCORDANCE {
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 16
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_CONCORDANCE"

    input:
    path ref_panel_with_index
    tuple val(sid), path(imputed_files), path(validate_files)
    tuple val(sid), path(concordance_txt)
    val ac_bins
    val min_val_gl 
    val min_val_dp

    output:
    tuple val(sid), path("${sid}.error.cal.txt.gz")  , emit: errors_cal
    tuple val(sid), path("${sid}.error.grp.txt.gz")  , emit: errors_grp
    tuple val(sid), path("${sid}.error.spl.txt.gz")  , emit: errors_spl
    tuple val(sid), path("${sid}.rsquare.grp.txt.gz"), emit: rsquare_grp
    tuple val(sid), path("${sid}.rsquare.spl.txt.gz"), emit: rsquare_spl
    tuple val(sid), path("${sid}_r2_sites.txt.gz")   , emit: rsquare_per_site, optional: true

    script:
    """
    GLIMPSE2_concordance            \
        --gt-val                    \
        --ac-bins $ac_bins          \
        --input $concordance_txt    \
        --thread $task.cpus         \
        --output $sid
    """
}
