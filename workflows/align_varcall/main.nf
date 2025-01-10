// Include processes
include { BWA_MEM                           } from '../../modules/bwa_mem'
include { SAMTOOLS_FLAGSTAT                 } from '../../modules/samtools/flagstat'
include { SAMTOOLS_INDEX                    } from '../../modules/samtools/index'
include { BCFTOOLS_MPILEUP                  } from '../../modules/bcftools/mpileup'
include { BCFTOOLS_INDEX                    } from '../../modules/bcftools/index'
include { BCFTOOLS_STATS as BCFTOOLS_STATS1 } from '../../modules/bcftools/stats'
include { MOSDEPTH                          } from '../../modules/mosdepth'

workflow ALIGN_VARCALL { 
    take:
    reference
    trimmed_reads
    bwaidx
    faidx

    main:
    BWA_MEM(trimmed_reads, reference, bwaidx)
    SAMTOOLS_FLAGSTAT(BWA_MEM.out.bam)
    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    MOSDEPTH(BWA_MEM.out.bam, reference)
    BCFTOOLS_MPILEUP(reference, BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai), faidx)
    BCFTOOLS_INDEX(BCFTOOLS_MPILEUP.out.bcf)
    BCFTOOLS_STATS1(BCFTOOLS_MPILEUP.out.bcf, 'before')

    emit:
    align = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
    bcfstats1 = BCFTOOLS_STATS1.out.bcfstats
}