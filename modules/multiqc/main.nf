// Define the `MULTIQC` process that performs report
process MULTIQC {
    container 'staphb/multiqc:latest'
    conda "${moduleDir}/environment.yml"
    tag "all_samples"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/MULTIQC"
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