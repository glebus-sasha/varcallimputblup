// Define the `COV_SUMMARY` process that merges all vcf files
process COV_SUMMARY {
    container ''
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/COV_SUMMARY"
    errorStrategy 'ignore'

    input:
    path(csvFiles)

    output:
    path("merged_${workflow.start.format('yyyyMMdd_HHmmss')}.csv"), emit: merged_csv

    script:
    """
    # Extract the header from the first file
    header=\$(head -n 1 "${csvFiles[0]}")

    # Create the merged file with a unique name
    merged_file="merged_${workflow.start.format('yyyyMMdd_HHmmss')}.csv"
    echo "\$header" > "\$merged_file"

    # Append the data from all files, skipping the header in subsequent files
    for file in "\${csvFiles[@]}"; do
        tail -n +2 "\$file" >> "\$merged_file"
    done
    """
}