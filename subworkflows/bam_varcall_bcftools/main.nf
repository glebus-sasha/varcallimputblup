// Include processes
include { BCFTOOLS_MPILEUP                          } from '../../modules/bcftools/mpileup'
include { BCFTOOLS_INDEX                            } from '../../modules/bcftools/index'
include { BCFTOOLS_STATS                            } from '../../modules/bcftools/stats'

workflow BAM_VARCALL_BCFTOOLS { 
    take:
    reference
    faidx
    align

    main:
    BCFTOOLS_MPILEUP(reference, align, faidx)
    BCFTOOLS_INDEX(BCFTOOLS_MPILEUP.out.bcf)
    BCFTOOLS_STATS(BCFTOOLS_MPILEUP.out.bcf, '')

    emit:
    bcf      = BCFTOOLS_MPILEUP.out.bcf
    csi      = BCFTOOLS_INDEX.out.csi
    bcfstats = BCFTOOLS_STATS.out.bcfstats
}