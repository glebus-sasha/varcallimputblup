// Define the `COV_SUMMARY` process that merges all vcf files
process COV_SUMMARY {
    container ''
    conda "${moduleDir}/environment.yml"
    tag ''
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_SUMMARY"
    errorStrategy 'ignore'

    input:
    path(csvFiles)

    output:
    path("merged_${workflow.start.format('yyyyMMdd_HHmmss')}.csv"), emit: merged_csv

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(readr)

    # Read all CSV files into a list of data frames
    csv_files <- c("${csvFiles.join('", "')}")
    data_frames <- lapply(csv_files, read_csv)

    # Combine all data frames into one
    combined_data <- bind_rows(data_frames)

    # Write the combined data to a new CSV file
    merged_file <- "merged_${workflow.start.format('yyyyMMdd_HHmmss')}.csv"
    write_csv(combined_data, merged_file)
    """
}