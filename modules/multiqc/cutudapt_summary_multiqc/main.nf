// Define the `CUTADAPT_SUMMARY_MULTIQC` process that performs report
process CUTADAPT_SUMMARY_MULTIQC {
    container 'staphb/multiqc:latest'
    conda "${moduleDir}/environment.yml"
    tag 'all samples'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/CUTADAPT_SUMMARY_MULTIQC"
//	  debug true
//    errorStrategy 'ignore'
	
    input:
    path fastqc1
    path fastqc2

    output:
    path '*.html', emit: html

    script:
    """
    multiqc .
    """
}