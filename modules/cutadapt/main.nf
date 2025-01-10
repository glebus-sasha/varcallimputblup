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
    tuple val(sid), path(read1), path(read2)
    tuple val(adapter5f), val(adapter5r)

    output:
    tuple val(sid), path("${sid}.R1.fq.gz"), path("${sid}.R2.fq.gz"), emit: cutadapted_reads
    
    script:
    """
    cutadapt -g $adapter5f -G $adapter5r -o ${sid}.R1.fq.gz -p ${sid}.R2.fq.gz ${read1} ${read2}
    """
}