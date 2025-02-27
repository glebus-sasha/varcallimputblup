include { BCF_TO_VCF                    } from '../../modules/bcftools/bcf_to_vcf'
include { VCF_TO_GDS                    } from '../../modules/local/vcf_to_gds'
include { COMBINE_GDS                   } from '../../modules/local/combine_gds'
include { GDS_CLUSTER                   } from '../../modules/local/gds_cluster'

workflow BCF_CLUSTERING{
    take:
    bcf

    main:
    BCF_TO_VCF(bcf)
    VCF_TO_GDS(BCF_TO_VCF.out.vcf)
    COMBINE_GDS(VCF_TO_GDS.out.gds.map{it -> it[1]}.collect())
    GDS_CLUSTER(COMBINE_GDS.out.combined_gds)

    emit:
    cluster_tab     = GDS_CLUSTER.out.cluster_tab
    dendrogram_nwk  = GDS_CLUSTER.out.dendrogram_nwk
    dendrogram_png  = GDS_CLUSTER.out.dendrogram_png
    pca_png         = GDS_CLUSTER.out.pca_png
}