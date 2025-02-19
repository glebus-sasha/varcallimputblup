// Define the `PREPARE_CONCORDANCE` process that plots concordance accuracy
process PREPARE_CONCORDANCE {
    conda ""
    tag { 
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
//  publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/PREPARE_CONCORDANCE"
//	debug true
    errorStrategy 'ignore'
	
    input:
    val chrs
    tuple val(sid), path(imputed_bcf), path(imputed_index)

    output:
    tuple val(sid), path("${sid}_concordance.txt"), emit: concordance_list
    tuple val(sid), path("*.bcf*"),                 emit: imputed

    script:
    def chr_list = chrs.collect { chr -> 
        """
            bcftools view -r ${chr} ${imputed_bcf} -o ${chr}_${imputed_bcf}
            bcftools index ${chr}_${imputed_bcf}
            echo "${chr} ${chr}.vcf.gz ${sid}_${chr}.bcf ${chr}_${imputed_bcf}" >> ${sid}_concordance.txt
        """
        }.join('\n')

    """
    ${chr_list}
    """
}