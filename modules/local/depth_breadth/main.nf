// Define the `DEPTH_BREADTH` process that calculates the depth of coverage and generates a plot and summary statistics
process DEPTH_BREADTH {
    container 'glebusasha/r_env_image:latest'
    conda "${moduleDir}/environment.yml"
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/DEPTH_BREADTH"
//	  debug true
    errorStrategy 'ignore'

    input:
    tuple val(sid), path(mosdepth_summary), path(coverage_width)
    path reference_length

    output:
    path "${sid}_stats.csv", emit: cov_stats

    script:
    """
    echo ${mosdepth_summary}
    echo ${reference_length}
    echo ${coverage_width}
    touch "${sid}_stats.csv"
    """
}
