// Include processes
include { FASTQC as FASTQC1                 } from '../../modules/fastqc'
include { FASTP                             } from '../../modules/fastp'
include { FASTQC as FASTQC2                 } from '../../modules/fastqc'
include { BWA_MEM                           } from '../../modules/bwa_mem'
include { SAMTOOLS_FLAGSTAT                 } from '../../modules/samtools/flagstat'
include { SAMTOOLS_INDEX                    } from '../../modules/samtools/index'
include { BCFTOOLS_MPILEUP                  } from '../../modules/bcftools/mpileup'
include { BCFTOOLS_INDEX                    } from '../../modules/bcftools/index'
include { BCFTOOLS_STATS as BCFTOOLS_STATS1 } from '../../modules/bcftools/stats'

workflow ALIGN_VARCALL { 
    take:
    reference
    input_fastqs
    bwaidx
    faidx

    main:
    FASTQC1(input_fastqs)
    FASTP(input_fastqs)
    FASTQC2(FASTP.out.trimmed_reads)
    BWA_MEM(FASTP.out.trimmed_reads, reference, bwaidx)
    SAMTOOLS_FLAGSTAT(BWA_MEM.out.bam)
    SAMTOOLS_INDEX(BWA_MEM.out.bam)
    BCFTOOLS_MPILEUP(reference, BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai), faidx)
    BCFTOOLS_INDEX(BCFTOOLS_MPILEUP.out.bcf)
    BCFTOOLS_STATS1(BCFTOOLS_MPILEUP.out.bcf, 'before')

    emit:
    fastp = FASTP.out.json
    align = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)
    fastqc1 = FASTQC1.out.zip
    fastqc2 = FASTQC2.out.zip
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
    bcfstats1 = BCFTOOLS_STATS1.out.bcfstats
}