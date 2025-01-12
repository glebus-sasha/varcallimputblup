// Define the `COV_STATS` process that calculates the depth of coverage and generates a plot and summary statistics
process COV_STATS {
    container 'glebusasha/r_env_image:latest'
    conda "${moduleDir}/environment.yml"
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_STATS"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(mosdepth_summary), path(coverage_width), path(bcfstatsFile)
    path reference_length

    output:
    tuple val(sid), path("${sid}_stats.csv")        , emit: cov_stats
    tuple val(sid), path("${sid}_coverage_plot.png"), emit: coverage_plot

    script:
    """
    #!/usr/bin/env Rscript
    source('${projectDir}/assets/process_chromosome_data.R')

    process_chromosome_data('$sid', '$mosdepth_summary', '$reference_length', '$coverage_width', '$bcfstatsFile')
    """
}
