// Define the `VARCALL_REF_BCFTOOLS` process that performs variant calling
process VARCALL_REF_BCFTOOLS {
    container 'staphb/bcftools:latest'
    conda 'bcftools'
    tag {
        sid.length() > 40 ? "${sid.take(20)}...${sid.takeRight(20)}" : sid
    }
    cpus 6
//  publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VARCALL_REF_BCFTOOLS"
//	debug true
//  errorStrategy 'ignore'
	
    input:
    path reference
    path fai
    tuple val(sid), path(bam), path(bai), val(chr), path(ref_panel), path(ref_panel_index), path(ref_panel_tsv), path(ref_panel_tsv_index)
    
    output:
    tuple val(sid), val(chr), path("${sid}_${chr}.bcf"), path("${sid}_${chr}.bcf.csi"), emit: bcf
    
    script:
    """    
    bcftools mpileup \
        -f $reference \
        -I -E -a 'FORMAT/DP' \
        -T ${ref_panel} \
        $bam \
        -Ou | \
        bcftools call \
        -Aim -C alleles -T ${ref_panel_tsv} \
        -Ob \
        -o ${sid}_${chr}.bcf \
        --threads ${task.cpus}
    bcftools index ${sid}_${chr}.bcf
    """
}
