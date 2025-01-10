// Define the `MOSDEPTH` process that calculates coverage depth
process MOSDEPTH {
    container 'biocontainers/mosdepth:0.3.2'
    conda 'mosdepth'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/MOSDEPTH"
//    debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bam)

    output:
    tuple val(sid), path("${sid}.mosdepth.global.dist.txt"), emit: global_dist
    tuple val(sid), path("${sid}.per-base.bed.gz")         , emit: per_base
    tuple val(sid), path("${sid}.regions.bed.gz")          , emit: regions

    script:
    """
        mosdepth -n -t 6 ${sid} ${bam}
    """
}
