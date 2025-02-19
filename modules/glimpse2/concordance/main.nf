process GLIMPSE2_CONCORDANCE {
    tag 'all_samples'
    cpus params.cpus
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_CONCORDANCE"

    input:
    path sids
    path ref_panel_with_index
    path imputed_bcf_files
    path merged_concordance
    path ch_validate_flat
    val ac_bins
    val min_val_gl 
    val min_val_dp

    output:
    tuple val(sid), path("*.error.cal.txt.gz")  , emit: errors_cal
    tuple val(sid), path("*.error.grp.txt.gz")  , emit: errors_grp
    tuple val(sid), path("*.error.spl.txt.gz")  , emit: errors_spl
    tuple val(sid), path("*.rsquare.grp.txt.gz"), emit: rsquare_grp
    tuple val(sid), path("*.rsquare.spl.txt.gz"), emit: rsquare_spl
    tuple val(sid), path("*_r2_sites.txt.gz")   , emit: rsquare_per_site, optional: true

    script:
    """
    GLIMPSE2_concordance            \
        --ac-bins $ac_bins          \
        --min-val-gl $min_val_gl    \
        --min-val-dp $min_val_dp    \
        --input $merged_concordance \
        --samples $sids             \
        --thread $task.cpus         \
        --output 'concordance'
    """
}
