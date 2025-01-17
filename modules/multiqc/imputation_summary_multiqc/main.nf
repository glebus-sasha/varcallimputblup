// Define the `IMPUTATION_SUMMARY_MULTIQC` process that performs report
process IMPUTATION_SUMMARY_MULTIQC {
    container 'staphb/multiqc:latest'
    conda "${moduleDir}/environment.yml"
    tag ""
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/IMPUTATION_SUMMARY_MULTIQC"
//	  debug true
//    errorStrategy 'ignore'
	
    input:
    path files

    output:
    path '*.html', emit: html

    script:
    """
    multiqc .
    """
}