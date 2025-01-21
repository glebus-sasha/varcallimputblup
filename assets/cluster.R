library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(dplyr)
library(stringr)
library(dendextend)
library(phylogram)
library(ggplot2)
library(readr)
library(ape)

# Read command line arguments
args <- commandArgs(TRUE)
input_dir <- args[1]
combined_gds_file <- args[2]
output_dir <- args[3]

# Function to convert VCF to GDS
convert_vcf_to_gds <- function(vcf_files, output_dir) {
  vcf_files %>%
    purrr::walk(~ {
      gds_file <- paste0(output_dir, "/", basename(.x), ".gds")
      snpgdsVCF2GDS(.x, gds_file, method = "biallelic.only")
    })
}

# Function to combine GDS files
combine_gds_files <- function(gds_files, combined_gds_file) {
  snpgdsCombineGeno(gds_files, combined_gds_file)
}

# Function to perform LD pruning and PCA
perform_pca <- function(gds_file, output_dir) {
  genofile <- snpgdsOpen(gds_file)
  
  # LD pruning with all chromosomes included
  set.seed(1000)
  snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, autosome.only=FALSE) # Adjust LD pruning threshold
  snpset.id <- unlist(snpset)
  
  # Extract sample IDs
  samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  samp.id <- str_remove_all(samp.id, ".*/|\\.sorted\\.bam")
  
  # Compute IBS matrix
  ibs_mat <- snpgdsIBS(genofile, num.thread = 2, autosome.only = FALSE)
  
  # Construct dendrogram
  ibs_hc <- snpgdsHCluster(ibs_mat)
  dend <- as.dendrogram(ibs_hc$dendrogram)
  
  # Save dendrogram
  dend_file <- paste0(output_dir, "/dendrogram.png")
  png(dend_file, width = 800, height = 1000)
  plot(dend)
  dev.off()
  
  # Perform PCA
  pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 40, autosome.only=FALSE)
  pc.percent <- pca$varprop * 100
  
  # Table with PCA results
  tab <- tibble(
    sample.id = samp.id,
    EV1 = pca$eigenvect[, 1],
    EV2 = pca$eigenvect[, 2]
  )
  
  # Save table
  write_tsv(tab, file = paste0(output_dir, "/cluster_tab.tsv"))
  
  # PCA visualization
  p <- ggplot(tab, aes(x = EV2, y = EV1)) +
    geom_point() +
    geom_text(aes(label = sample.id)) +
    theme_minimal()
  
  # Save plot
  ggsave(paste0(output_dir, "/cluster.png"), plot = p, width = 15, height = 6, units = "in", dpi = 300)
  
  # Close GDS file
  snpgdsClose(genofile)
}

# Main function to execute all steps
main <- function(input_dir, combined_gds_file, output_dir) {
  vcf_files <- list.files(path = input_dir, pattern = "\\.vcf$", full.names = TRUE)
  convert_vcf_to_gds(vcf_files, output_dir)
  
  gds_files <- list.files(path = output_dir, pattern = "\\.gds$", full.names = TRUE)
  combine_gds_files(gds_files, combined_gds_file)
  
  perform_pca(combined_gds_file, output_dir)
}

# Execute main function
main(input_dir, combined_gds_file, output_dir)