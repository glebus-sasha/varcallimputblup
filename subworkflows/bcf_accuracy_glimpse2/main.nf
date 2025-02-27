include { BCFTOOLS_INDEX                        } from '../../modules/bcftools/index'
include { BCFTOOLS_STATS                        } from '../../modules/bcftools/stats'
include { GLIMPSE2_CONCORDANCE                  } from '../../modules/glimpse2/concordance'
include { BCF_RENAME_BCFTOOLS                   } from '../../modules/local/bcf_rename_bcftools'
include { PREPARE_CONCORDANCE                   } from '../../modules/local/prepare_concordance'
include { MERGE_PREPARE_CONCORDANCE             } from '../../modules/local/merge_prepare_concordance'
include { VARCALL_REF_BCFTOOLS                  } from '../../modules/local/varcall_ref_bcftools'
include { IMPUTATION_ACCURACY_PLOT              } from '../../modules/local/imputation_accuracy_plot'

workflow BCF_ACCURACY_GLIMPSE2{
    take:
    reference
    faidx
    imputed_bcf
    align
    ref_panel

    main:
    BCFTOOLS_INDEX(imputed_bcf)
    VARCALL_REF_BCFTOOLS(reference, faidx, align.combine(ref_panel))

    ch_validate = VARCALL_REF_BCFTOOLS.out.bcf
    ch_imputed  = imputed_bcf.join(BCFTOOLS_INDEX.out.csi)
    
    BCF_RENAME_BCFTOOLS(ch_validate.join(ch_imputed))
   
    ch_chrs     = ref_panel.map{it[0]}.collect()
   
    PREPARE_CONCORDANCE(ch_chrs, BCF_RENAME_BCFTOOLS.out.imputed_fixed_bcf)
    
    ch_ref_panel_with_index = ref_panel.map{[it[1], it[2]]}.collect()
    ch_i                    = PREPARE_CONCORDANCE.out.imputed
    ch_v                    = ch_validate.map { [it[0], [it[2], it[3]]] }
        .groupTuple()
        .map { sid, files -> [sid, files.flatten()] }
    ch_i_v                  = ch_i.join(ch_v)
    ch_concordance_txt      = PREPARE_CONCORDANCE.out.concordance_list

    GLIMPSE2_CONCORDANCE(
        ch_ref_panel_with_index, 
        ch_i_v, 
        ch_concordance_txt, 
        '1 5 10 20 50 100 200 500 1000 2000 5000 10000  20000 50000 100000 15000', 
        0, 0
        )

    IMPUTATION_ACCURACY_PLOT(GLIMPSE2_CONCORDANCE.out.rsquare_grp)

    emit:
    glimpse2_errors_grp = GLIMPSE2_CONCORDANCE.out.errors_grp
    glimpse2_errors_spl = GLIMPSE2_CONCORDANCE.out.errors_spl
}