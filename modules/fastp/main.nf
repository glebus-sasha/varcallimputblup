// Define the `FASTP` process that performs quality trimming and filtering of reads
process FASTP{
    
    container = 'nanozoo/fastp:0.23.1--9f2e255'
    tag "${sid}"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/FASTP"
    cpus 10
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(read1), path(read2)

    output:
    tuple val(sid), path("${sid}.R1.fastq.gz"), path("${sid}.R2.fastq.gz"), emit: trimmed_reads
    path '*.html', emit: html, optional: true
    path '*.json', emit: json, optional: true

    script:
    """
    fastp \
    --thread ${task.cpus} \
    --in1 $read1 \
    --in2 $read2\
    --out1 "${sid}.R1.fastq.gz" \
    --out2 "${sid}.R2.fastq.gz" \
    --html ${sid}.fastp_stats.html \
    --json ${sid}.fastp_stats.json 
    """
}
