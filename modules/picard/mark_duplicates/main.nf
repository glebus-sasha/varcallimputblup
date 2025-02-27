process PICARD_MARK_DUPLICATES {
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    conda 'bioconda::picard=2.27.4'
    container 'broadinstitute/picard:latest'
    //publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/PICARD_MARK_DUPLICATES"


    input:
    tuple val(sid), path(bam), path(bai)

    output:
    tuple val(sid), path("${sid}.dedup.bam")        , emit: bam
    tuple val(sid), path("${sid}.dedup.metrics.txt"), emit: metrics

    script:
    """
    picard MarkDuplicates \
        I=$bam \
        O=${sid}.dedup.bam \
        M=${sid}.dedup.metrics.txt \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=LENIENT
    """
}