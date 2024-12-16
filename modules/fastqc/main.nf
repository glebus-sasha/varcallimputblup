// Define the `FASTQC` process that performs quality trimming and filtering of reads
process FASTQC {
    container = 'staphb/fastqc:0.12.1'
    tag "${sid}"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/FASTQC", pattern: '*.html'
//	debug true
//  errorStrategy 'ignore'

    input:
    tuple val(sid), path(read1), path(read2)
    val tag

    output:
    path "*.html", emit: html
    path "*.zip" , emit: zip

    script:
    """
    fastqc $read1 $read2 --threads 6
    mv ${read1}_fastqc.html > ${read1}_${tag}_fastqc.html
    mv ${read2}_fastqc.html > ${read1}_${tag}_fastqc.html

    """

    stub:
    """

    """
}
