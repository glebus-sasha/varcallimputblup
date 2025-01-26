// Include processes
include { GLIMPSE2_CHUNK                     } from '../../modules/glimpse2/chunk'
include { GLIMPSE2_CONCORDANCE               } from '../../modules/glimpse2/concordance'
include { GLIMPSE2_LIGATE                    } from '../../modules/glimpse2/ligate'
include { GLIMPSE2_PHASE                     } from '../../modules/glimpse2/phase'
include { GLIMPSE2_SPLITREFERENCE            } from '../../modules/glimpse2/splitreference'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_2 } from '../../modules/bcftools/stats'

workflow IMPUTE {
    take:
    ref_panel_with_index
    ref_panel_index
    align

    main:
    GLIMPSE2_CHUNK(ref_panel_with_index)
    chr_IRG_ORG = GLIMPSE2_CHUNK.out.chunk_chr.flatMap { chr, file ->
        file.splitCsv(header:false,sep:'\t').collect { coord ->
            [chr, coord[2], coord[3]]
        }
    }
    GLIMPSE2_SPLITREFERENCE(ref_panel_with_index.combine(chr_IRG_ORG, by: 0))
    GLIMPSE2_PHASE(
        GLIMPSE2_SPLITREFERENCE.out.bin_ref.combine(align).combine(ref_panel_index, by: 0)
    )
    GLIMPSE2_LIGATE(GLIMPSE2_PHASE.out.phased_variants.groupTuple())
    BCFTOOLS_STATS_2(GLIMPSE2_LIGATE.out.merged_variants, 'after')

    emit:
    imputed_bcf = GLIMPSE2_LIGATE.out.merged_variants
    bcfstats_imputed = BCFTOOLS_STATS_2.out.bcfstats
}