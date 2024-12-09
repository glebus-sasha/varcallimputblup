// Define the `FASTP` process that performs quality trimming and filtering of reads
process FASTP{
    
    container = 'nanozoo/fastp:0.23.1--9f2e255'
    tag "${sid}"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/FASTP"
    cpus 10
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(reads)

    output:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), emit: trimmed_reads
    path '*.html', emit: html, optional: true
    path '*.json', emit: json, optional: true

    script:
    fq_1_trimmed = sid + '_R1.fastq.gz'
    fq_2_trimmed = sid + '_R2.fastq.gz'
    """
    fastp \
    --thread ${task.cpus} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]}\
    --out1 $fq_1_trimmed \
    --out2 $fq_2_trimmed \
    --html ${sid}.fastp_stats.html \
    --json ${sid}.fastp_stats.json 
    """
}
