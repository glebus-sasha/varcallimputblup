// Define the `COV_SUMMARY` process that generates a summary table for all samples
process COV_SUMMARY {
    container ''
    conda "${moduleDir}/environment.yml"
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 10
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_SUMMARY"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(breadthFile), path(depthStatsFile), path(bcfstatsFile)

    output:
    path("${sid}_stats.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    # Load the R script
    source('${projectDir}/assets/coverage_summary.R')

    # Call the function with the provided inputs
    create_and_save_summary_table('${sid}', '${depthStatsFile}', '${bcfstatsFile}', '${breadthFile}')
    """
}