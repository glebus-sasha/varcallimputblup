// Include processes
include { GLIMPSE2_CHUNK          } from '../../modules/glimpse2/chunk'
include { GLIMPSE2_CONCORDANCE    } from '../../modules/glimpse2/concordance'
include { GLIMPSE2_LIGATE         } from '../../modules/glimpse2/ligate'
include { GLIMPSE2_PHASE          } from '../../modules/glimpse2/phase'
include { GLIMPSE2_SPLITREFERENCE } from '../../modules/glimpse2/splitreference'
include { BCFTOOLS_INDEX          } from '../../modules/bcftools/index'
include { BCFTOOLS_STATS          } from '../../modules/bcftools/stats'

workflow BAM_IMPUTE_GLIMPSE2 {
    take:
    ref_panel
    align

    main:
    GLIMPSE2_CHUNK(ref_panel)

    chr_IRG_ORG = GLIMPSE2_CHUNK.out.chunk_chr.flatMap { chr, file ->
        file.splitCsv(header:false,sep:'\t').collect { coord ->
            [chr, coord[2], coord[3]]
        }
    }

    GLIMPSE2_SPLITREFERENCE(ref_panel.combine(chr_IRG_ORG, by: 0))
    GLIMPSE2_PHASE(
        GLIMPSE2_SPLITREFERENCE.out.bin_ref.combine(align)
    )
    GLIMPSE2_LIGATE(GLIMPSE2_PHASE.out.phased_variants.groupTuple())
    BCFTOOLS_INDEX(GLIMPSE2_LIGATE.out.merged_variants)
    BCFTOOLS_STATS(GLIMPSE2_LIGATE.out.merged_variants.join(BCFTOOLS_INDEX.out.csi), '_imputed')

    emit:
    bcf              = GLIMPSE2_LIGATE.out.merged_variants
    csi              = BCFTOOLS_INDEX.out.csi
    bcfstats_imputed = BCFTOOLS_STATS.out.bcfstats
}