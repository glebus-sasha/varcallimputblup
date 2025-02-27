// Include processes
include { BCFTOOLS_VIEW } from '../../modules/bcftools/view'

workflow CARPET { 
    take:
    bcf_csi
    bed

    main:
    BCFTOOLS_VIEW(bcf_csi, bed, "")


    emit:
    bcf = ""
}