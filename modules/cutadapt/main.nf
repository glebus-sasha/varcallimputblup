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
    tuple val(adapter_3f), val(adapter_3r), val(adapter5f), val(adapter5r)

    output:
    tuple val(sid), path("${sid}.R1.fq.gz"), path("${sid}.R2.fq.gz"), emit: cutadapted_reads
    

    script:
    def command = "cutadapt"

    if (adapter_3f) {
        command += " -a ${adapter_3f}"
    }
    if (adapter_3r) {
        command += " -A ${adapter_3r}"
    }
    if (adapter5f) {
        command += " -g ${adapter5f}"
    }
    if (adapter5r) {
        command += " -G ${adapter5r}"
    }

    command += " -o ${sid}.R1.fq.gz -p ${sid}.R2.fq.gz ${read1} ${read2}"
    """
    ${command}
    """
}