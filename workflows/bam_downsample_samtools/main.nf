// Include processes
include { SAMTOOLS_VIEW as SAMTOOLS_DOWNSAMPLE } from '../../modules/samtools/view'
include { SAMTOOLS_INDEX                       } from '../../modules/samtools/index'

workflow BAM_DOWNSAMPLE_SAMTOOLS { 
    take:
    bam

    main:
    SAMTOOLS_DOWNSAMPLE(bam, 0.9)
    SAMTOOLS_INDEX(SAMTOOLS_DOWNSAMPLE.out.bam)

    emit:
    align = SAMTOOLS_DOWNSAMPLE.out.bam.join(SAMTOOLS_INDEX.out.bai)
}