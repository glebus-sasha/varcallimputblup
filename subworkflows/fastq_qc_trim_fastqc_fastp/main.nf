// Include processes
include { FASTQC              } from '../../modules/fastqc'
include { FASTP               } from '../../modules/fastp'
include { FASTQC as FASTQC_2  } from '../../modules/fastqc'

workflow FASTQ_QC_TRIM_FASTQ_FASTP { 
    take:
    input_fastqs

    main:
    FASTQC(input_fastqs)
    FASTP(input_fastqs)
    FASTQC_2(FASTP.out.fastq_gz)

    emit:
    trimmed_reads   = FASTP.out.fastq_gz
    fastp_json      = FASTP.out.json
    fastqc          = FASTQC.out.zip
    fastqc_trimmed  = FASTQC_2.out.zip
}