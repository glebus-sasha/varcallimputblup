// Include processes
include { FASTQC as FASTQC_BEFORE } from '../../modules/fastqc'
include { FASTP                   } from '../../modules/fastp'
include { CUTADAPT                } from '../../modules/cutadapt'
include { FASTQC as FASTQC_AFTER  } from '../../modules/fastqc'

workflow CUTADAPT_QC{
    take:
    input_fastqs
    adapters
    
    main:
    FASTQC_BEFORE(input_fastqs)
    CUTADAPT(input_fastqs, adapters)
    FASTQC_AFTER(FASTP.out.trimmed_reads)    
}