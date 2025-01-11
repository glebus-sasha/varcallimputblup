// Define the `GENOME_LENGTH` process that calculates the length of a reference sequence from a FASTA file
process GENOME_LENGTH {
    container ''
    conda ""
    tag "$reference"
    // publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GENOME_LENGTH"
    // debug true
    // errorStrategy 'ignore'

    input:
    path reference

    output:
    path "${reference.baseName}_length.txt", emit: genome_length

    script:
    """
    `grep -v '^>' "$reference" | tr -d '\n' | wc -c` > "${reference.baseName}_length.txt"
    """
}