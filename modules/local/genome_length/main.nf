// Define the `GENOME_LENGTH` process that calculates the length of a reference sequence from a FASTA file
process GENOME_LENGTH {
    container ''
    conda ""
    tag {
        reference.length() > 40 ? "${reference.take(20)}...${reference.takeRight(20)}" : reference
    }
    // publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/GENOME_LENGTH"
    // debug true
    errorStrategy 'ignore'

    input:
    path reference

    output:
    val length, emit: genome_length

    script:
    """
    length=$(grep -v '^>' "$reference" | tr -d '\n' | wc -c)
    """
}