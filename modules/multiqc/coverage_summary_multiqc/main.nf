// Define the `COVERAGE_SUMMARY_MULTIQC` process that performs report
process COVERAGE_SUMMARY_MULTIQC {
    container 'staphb/multiqc:latest'
    conda "${moduleDir}/environment.yml"
    tag 'all samples'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}"
//	  debug true
//    errorStrategy 'ignore'
	
    input:
    path fastp
    path fastqc1
    path fastqc2
    path flagstat
    path bcfstats1
    path mosdepth

    output:
    path '*.html', emit: html

    script:
    """
    multiqc .
    """
}