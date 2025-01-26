// Include processes
include { BWA_MEM                                   } from '../../modules/bwa_mem'
include { SAMTOOLS_RMDUP                            } from '../../modules/samtools/rmdup'
include { SAMTOOLS_FLAGSTAT                         } from '../../modules/samtools/flagstat'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_2  } from '../../modules/samtools/flagstat'
include { SAMTOOLS_INDEX                            } from '../../modules/samtools/index'
include { BCFTOOLS_MPILEUP                          } from '../../modules/bcftools/mpileup'
include { BCFTOOLS_INDEX                            } from '../../modules/bcftools/index'
include { BCFTOOLS_STATS                            } from '../../modules/bcftools/stats'
include { MOSDEPTH                                  } from '../../modules/mosdepth'

workflow ALIGN_VARCALL { 
    take:
    reference
    trimmed_reads
    bwaidx
    faidx

    main:
    BWA_MEM(trimmed_reads, reference, bwaidx)
    SAMTOOLS_RMDUP(BWA_MEM.out.bam)
    SAMTOOLS_FLAGSTAT(BWA_MEM.out.bam, '')
    SAMTOOLS_INDEX(SAMTOOLS_RMDUP.out.bam)
    SAMTOOLS_FLAGSTAT_2(SAMTOOLS_RMDUP.out.bam, '_dedup')

    MOSDEPTH(
        SAMTOOLS_RMDUP.out.bam          |
        join(SAMTOOLS_INDEX.out.bai)
    )

    BCFTOOLS_MPILEUP(
        reference, 
        SAMTOOLS_RMDUP.out.bam.join(SAMTOOLS_INDEX.out.bai), 
        faidx
        )
    BCFTOOLS_INDEX(BCFTOOLS_MPILEUP.out.bcf)
    BCFTOOLS_STATS(BCFTOOLS_MPILEUP.out.bcf, '')

    emit:
    bcf                 = BCFTOOLS_MPILEUP.out.bcf
    align               = SAMTOOLS_RMDUP.out.bam.join(SAMTOOLS_INDEX.out.bai)
    flagstat            = SAMTOOLS_FLAGSTAT.out.flagstat
    flagstat_dedup      = SAMTOOLS_FLAGSTAT_2.out.flagstat
    bcfstats            = BCFTOOLS_STATS.out.bcfstats
    mosdepth            = MOSDEPTH.out.global_dist
    mosdepth_summary    = MOSDEPTH.out.summary
}