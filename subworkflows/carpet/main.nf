// Include processes
include { BCFTOOLS_VIEW                     } from '../../modules/bcftools/view'
include { VCF_CARPET                        } from '../../modules/local/vcf_carpet'
include { VCF_CARPET as VCF_CARPET_1        } from '../../modules/local/vcf_carpet'
include { BCF_TO_VCF                        } from '../../modules/bcftools/bcf_to_vcf'
include { BCFTOOLS_MPILEUP                  } from '../../modules/bcftools/mpileup'
include { BCFTOOLS_INDEX                    } from '../../modules/bcftools/index'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_1  } from '../../modules/bcftools/view'

workflow CARPET { 
    take:
    bcf_csi
    bed
    interpretation
    align
    reference
    faidx

    main:
    BCFTOOLS_VIEW(bcf_csi, bed, "_imputed")
    /*VCF_CARPET(
        BCFTOOLS_VIEW.out.vcf.map{it[1]}.collect(), 
        bed, 
        interpretation,
        'imputed'
        )*/
    BCFTOOLS_MPILEUP(reference, align, faidx)
    BCFTOOLS_INDEX(BCFTOOLS_MPILEUP.out.bcf)
    BCFTOOLS_VIEW_1(BCFTOOLS_MPILEUP.out.bcf.join(BCFTOOLS_INDEX.out.csi), bed, '_original')
    //BCF_TO_VCF(BCFTOOLS_MPILEUP.out.bcf.join(BCFTOOLS_INDEX.out.csi))
    VCF_CARPET_1(
        BCFTOOLS_VIEW_1.out.vcf.map{it[1]}.mix(BCFTOOLS_VIEW.out.vcf.map{it[1]}).collect(), 
        bed, 
        interpretation,
        'original'
        )

    emit:
    vcf         = BCFTOOLS_VIEW.out.vcf
    /*carpet      = VCF_CARPET.out.carpet
    table_csv   = VCF_CARPET.out.table_csv
    table_html  = VCF_CARPET.out.table_html*/
}