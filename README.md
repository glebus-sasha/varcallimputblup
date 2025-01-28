# varcallimputblup
Nextflow pipeline for variant calling, imputation and cattle model construction

```
nextflow run -latest glebus-sasha/varcallimputblup -params-file params.yaml -profile mamba
```

```mermaid
%%{init: {'theme':'base'}}%%
flowchart TB
    subgraph " "
    v0["Channel.fromPath"]
    v2["Channel.fromFilePairs"]
    v3["Channel.fromPath"]
    v5["Channel.fromPath"]
    v7["Channel.fromPath"]
    v9["Channel.fromPath"]
    v20["tag"]
    v24["tag"]
    v29["tag"]
    v32["tag"]
    v38["tag"]
    v73["report_title"]
    end
    subgraph " "
    v12["ref_panel_with_index"]
    v14[" "]
    v16[" "]
    v18[" "]
    v26["mosdepth_summary"]
    v37[" "]
    v49[" "]
    v56["pca_png"]
    v57["dendrogram_png"]
    v58["dendrogram_nwk"]
    v59["cluster_tab"]
    v75[" "]
    end
    subgraph clustering
    subgraph CLUSTERING
    subgraph QC_TRIM
    v13([FASTQC])
    v15([FASTP])
    v17([FASTQC_2])
    end
    subgraph ALIGN_VARCALL
    v19([BWA_MEM])
    v21([SAMTOOLS_FLAGSTAT])
    v22([SAMTOOLS_INDEX_2])
    v25([MOSDEPTH])
    v27([SAMTOOLS_RMDUP])
    v28([SAMTOOLS_INDEX])
    v30([SAMTOOLS_FLAGSTAT_2])
    v33([MOSDEPTH_2])
    v35([BCFTOOLS_MPILEUP])
    v36([BCFTOOLS_INDEX])
    v39([BCFTOOLS_STATS])
    v1(( ))
    v4(( ))
    v6(( ))
    v23(( ))
    v31(( ))
    v34(( ))
    end
    subgraph COVERAGE_SUMMARY
    v41([BAM_BREADTH])
    v42([GENOME_LENGTH])
    v45([COV_STATS])
    v48([COV_SUMMARY])
    v40(( ))
    v43(( ))
    v46(( ))
    end
    subgraph BCF_CLUSTERING
    v50([BCF_TO_VCF])
    v51([VCF_TO_GDS])
    v54([COMBINE_GDS])
    v55([GDS_CLUSTER])
    v52(( ))
    end
    v74([MULTIQC])
    v60(( ))
    end
    end
    v8(( ))
    v0 --> v1
    v2 --> v13
    v2 --> v15
    v3 --> v4
    v5 --> v6
    v7 --> v8
    v9 --> v8
    v8 --> v12
    v13 --> v14
    v13 --> v60
    v15 --> v17
    v15 --> v16
    v15 --> v19
    v15 --> v60
    v17 --> v18
    v17 --> v60
    v1 --> v19
    v4 --> v19
    v19 --> v21
    v19 --> v22
    v19 --> v27
    v19 --> v23
    v20 --> v21
    v21 --> v60
    v22 --> v23
    v24 --> v25
    v23 --> v25
    v25 --> v26
    v25 --> v60
    v27 --> v28
    v27 --> v30
    v27 --> v31
    v27 --> v34
    v27 --> v40
    v28 --> v31
    v28 --> v34
    v28 --> v40
    v29 --> v30
    v30 --> v60
    v32 --> v33
    v31 --> v33
    v33 --> v43
    v33 --> v60
    v1 --> v35
    v6 --> v35
    v34 --> v35
    v35 --> v36
    v35 --> v39
    v35 --> v50
    v36 --> v37
    v38 --> v39
    v39 --> v43
    v39 --> v60
    v40 --> v41
    v41 --> v43
    v1 --> v42
    v42 --> v45
    v43 --> v45
    v45 --> v46
    v46 --> v48
    v48 --> v49
    v50 --> v51
    v51 --> v52
    v52 --> v54
    v54 --> v55
    v55 --> v59
    v55 --> v58
    v55 --> v57
    v55 --> v56
    v73 --> v74
    v60 --> v74
    v74 --> v75
```

FASTQC: Quality control of raw sequencing data using FastQC.

FASTP: Trimming of reads to remove adapters and low-quality sequences using fastp.

BWA_MEM: Alignment of reads to the reference genome using BWA MEM.

SAMTOOLS_FLAGSTAT: Quality assessment of alignment using SAMtools flagstat.

SAMTOOLS_INDEX: Indexing of BAM files using SAMtools index.

MOSDEPTH: Calculation of sequencing depth using mosdepth.

SAMTOOLS_RMDUP: Removal of PCR duplicates using SAMtools rmdup.

BCFTOOLS_MPILEUP: Variant calling using BCFtools mpileup.

BCFTOOLS_INDEX: Indexing of VCF files using BCFtools index.

BCFTOOLS_STATS: Statistical analysis of variant calls using BCFtools stats.

BAM_BREADTH: Calculation of breadth of coverage using BAM files.

GENOME_LENGTH: Calculation of genome length.

COV_STATS: Statistical analysis of coverage using coverage statistics.

COV_SUMMARY: Summary of coverage statistics.

BCF_TO_VCF: Conversion of BCF files to VCF format.

VCF_TO_GDS: Conversion of VCF files to GDS format.

COMBINE_GDS: Combining multiple GDS files.

GDS_CLUSTER: Clustering of GDS files.
MULTIQC: Compilation of a comprehensive report including QC metrics, alignment results, and variant calling statistics.
