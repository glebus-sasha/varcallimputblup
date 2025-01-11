// Define the `DEPTH_BREADTH` process that calculates the depth of coverage and generates a plot and summary statistics
process BAM_DEPTH {
    container 'glebusasha/r_env_image:latest'
    conda "${moduleDir}/environment.yml"
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/DEPTH_BREADTH"
//	  debug true
    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile), path(bamIndex)

    output:
    tuple val(sid), path("${sid}_depth_plot.png") , emit: depth_plot
    tuple val(sid), path("${sid}_depth_stats.csv"), emit: depth_stats

    script:
    """

    """
}
