// Define the `BCF_RENAME_BCFTOOLS` process that renames sampl names in bcf
process BCF_RENAME_BCFTOOLS {
    container ''
    conda 'bcftools'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 10
//    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BCF_RENAME_BCFTOOLS"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), val(chr), path(validate_bcf), path(validate_index), path(imputed_bcf), path(imputed_index)

    output:
    tuple val(sid), path("${sid}_imputed_fixed.bcf"), path("${sid}_imputed_fixed.bcf.csi"), emit: imputed_fixed_bcf

    script:
    """
    TRUTH_SAMPLE=\$(bcftools query -l "$validate_bcf")
    IMPUTED_SAMPLE=\$(bcftools query -l "$imputed_bcf")

    echo -e "\$IMPUTED_SAMPLE\t\$TRUTH_SAMPLE" > rename_samples.txt
    bcftools reheader -s rename_samples.txt -o "${sid}_imputed_fixed.bcf" "$imputed_bcf" 
    bcftools index -f "${sid}_imputed_fixed.bcf"
    """
}