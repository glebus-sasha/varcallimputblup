// Define the `VCF_CLUSTER` process that performs clustering
process VCF_CLUSTER {
    conda "${moduleDir}/environment.yml"
    tag 'all_samples'
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VCF_CLUSTER"
//	debug true
    errorStrategy 'ignore'
	
    input:
    path bcf_files

    output:
    path 'cluster.png'
    path 'cluster_tab.tsv'
    path 'dendrogram.png'
    path 'dendrogram.nwk'

    script:
    """    
    mkdir .temp
    # Перебираем все файлы в переменной $bcf_files
    for bcf_file in $bcf_files; do
        # Определяем имя выходного файла .vcf
        vcf_file=".temp/\$(basename \$bcf_file .bcf).vcf"

        # Преобразуем .bcf в .vcf с помощью bcftools
        bcftools view \$bcf_file -O v -o \$vcf_file
    done
    Rscript ${projectDir}/assets/cluster.R ./.temp combined.gds ./
    """
}