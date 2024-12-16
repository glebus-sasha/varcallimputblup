// Include processes
include { GLIMPSE2_CHUNK                    } from '../../modules/glimpse2/chunk'
include { GLIMPSE2_CONCORDANCE              } from '../../modules/glimpse2/concordance'
include { GLIMPSE2_LIGATE                   } from '../../modules/glimpse2/ligate'
include { GLIMPSE2_PHASE                    } from '../../modules/glimpse2/phase'
include { GLIMPSE2_SPLITREFERENCE           } from '../../modules/glimpse2/splitreference'
include { BCFTOOLS_STATS as BCFTOOLS_STATS2 } from '../../modules/bcftools/stats'
include { MULTIQC                           } from '../../modules/multiqc'

workflow IMPUTE {
    take:
    ref_panel_with_index
    ref_panel_index
    align

    main:
    GLIMPSE2_CHUNK(ref_panel_with_index)
    IRG_ORG = GLIMPSE2_CHUNK.out.chunk_chr.splitCsv(header:false,sep:'\t').map{coord->[coord[2],coord[3]]}
    GLIMPSE2_SPLITREFERENCE(ref_panel_with_index.combine(IRG_ORG))
    GLIMPSE2_PHASE(
        align.combine(GLIMPSE2_SPLITREFERENCE.out.bin_ref.map{it->it[1]}).combine(ref_panel_index)
        )
    GLIMPSE2_LIGATE(GLIMPSE2_PHASE.out.phased_variants.groupTuple())
    BCFTOOLS_STATS2(GLIMPSE2_LIGATE.out.merged_variants, 'after')

    emit:
    imputed_bcf = GLIMPSE2_LIGATE.out.merged_variants
    bcfstats2 = BCFTOOLS_STATS2.out.bcfstats
}