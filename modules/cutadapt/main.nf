// Define the `CUTADAPT` process that performs variant calling
process CUTADAPT {
    container 'zavolab/cutadapt:1.16-slim'
    conda 'cutadapt'
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/CUTADAPT"
//	  debug true
//    errorStrategy 'ignore'
	
    input:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed)
    val adapter

    output:
    tuple val(sid), path("${sid}_R1.fq.gz"), path("${sid}_R2.fq.gz"), emit: cutadapted_reads
    
    script:
    """
    cutadapt -g $adapter -G $adapter -o ${sid}_R1.fq.gz -p ${sid}_R2.fq.gz ${fq_1_trimmed} ${fq_2_trimmed}
    """
}