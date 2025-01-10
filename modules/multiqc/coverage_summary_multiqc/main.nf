// Define the `COVERAGE_SUMMARY_MULTIQC` process that performs report
process COVERAGE_SUMMARY_MULTIQC {
    container 'staphb/multiqc:latest'
    conda 'multiqc python=3.12'
    tag 'all samples'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COVERAGE_SUMMARY_MULTIQC"
//	  debug true
//    errorStrategy 'ignore'
	
    input:
    path fastp
    path fastqc1
    path fastqc2
    path flagstat
    path bcfstats1

    output:
    path '*.html', emit: html

    script:
    """
    multiqc .
    """
}