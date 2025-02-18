// Include processes
include { BWA_MEM                                   } from '../../modules/bwa_mem'
include { SAMTOOLS_FLAGSTAT                         } from '../../modules/samtools/flagstat'
include { SAMTOOLS_INDEX                            } from '../../modules/samtools/index'
include { MOSDEPTH                                  } from '../../modules/mosdepth'

workflow FASTQ_ALIGN_BWA { 
    take:
    reference
    fastq
    bwaidx

    main:
    BWA_MEM(fastq, reference, bwaidx)
    SAMTOOLS_FLAGSTAT(BWA_MEM.out.bam, '')
    SAMTOOLS_INDEX(BWA_MEM.out.bam)

    align = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)
    MOSDEPTH(align, '')

    emit:
    align            = align
    flagstat         = SAMTOOLS_FLAGSTAT.out.flagstat
    mosdepth         = MOSDEPTH.out.global_dist
    mosdepth_summary = MOSDEPTH.out.summary
}