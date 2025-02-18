process GLIMPSE2_CONCORDANCE {
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container 'imary116/glimpse2:with-bcftools-and-updated-info-score'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GLIMPSE2_CONCORDANCE"

    input:
    val chrs
    path ref_panel_with_index
    tuple val(sid), path(imputed_bcf), path(imputed_index), path(validate_files) 

    output:
    tuple val(sid), path("*.error.cal.txt.gz")  , emit: errors_cal
    tuple val(sid), path("*.error.grp.txt.gz")  , emit: errors_grp
    tuple val(sid), path("*.error.spl.txt.gz")  , emit: errors_spl
    tuple val(sid), path("*.rsquare.grp.txt.gz"), emit: rsquare_grp
    tuple val(sid), path("*.rsquare.spl.txt.gz"), emit: rsquare_spl
    tuple val(sid), path("*_r2_sites.txt.gz")   , emit: rsquare_per_site, optional: true

    script:
    def chr_list = chrs.collect { chr -> 
        """
        bcftools view -r ${chr} ${imputed_bcf} -o ${chr}_${imputed_bcf}
        bcftools index ${chr}_${imputed_bcf}
        echo "${chr} ${chr}.vcf.gz ${sid}_${chr}.bcf ${chr}_${imputed_bcf}" >> concordance.txt
        """
    }.join('\n')

    """
    ${chr_list}

    GLIMPSE2_concordance \\
        --ac-bins 1 5 10 20 \\
        --min-val-gl 0.1 \\
        --min-val-dp 1 \\
        --input concordance.txt \\
        --thread $task.cpus \\
        --output $sid
    """
}
