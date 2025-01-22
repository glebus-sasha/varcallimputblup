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
library(ggfortify)

# Function to perform LD pruning, PCA, and return both dendrogram and PCA results
gds_cluster <- function() {
  # Automatically find the combined.gds file in the current directory
  gds_file <- "combined.gds"
  if (!file.exists(gds_file)) {
    stop("combined.gds file not found in the current directory.")
  }

  genofile <- snpgdsOpen(gds_file)

  # LD pruning with all chromosomes included
  set.seed(1000)
  snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, autosome.only=FALSE)
  snpset.id <- unlist(snpset)

  # Extract sample IDs
  samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
  samp.id <- str_remove_all(samp.id, ".*/|\\.sorted\\.bam")

  # Compute IBS matrix
  ibs_mat <- snpgdsIBS(genofile, num.thread = 2, autosome.only = FALSE)

  # Construct dendrogram
  ibs_hc <- snpgdsHCluster(ibs_mat)
  dend <- as.dendrogram(ibs_hc$dendrogram)

  # Save dendrogram as Newick format
  newick_file <- "dendrogram.nwk"
  write.tree(as.phylo(dend), file = newick_file)

  # Save dendrogram
  dend_file <- "dendrogram.png"
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
  write_tsv(tab, file = "cluster_tab.tsv")

  # PCA visualization
  p <- ggplot(tab, aes(x = EV2, y = EV1)) +
    geom_point() +
    geom_text(aes(label = sample.id)) +
    theme_minimal()

  # Save PCA plot
  pca_file <- "pca.png"
  ggsave(pca_file, plot = p, width = 15, height = 6, units = "in", dpi = 300)

  # Close GDS file
  snpgdsClose(genofile)

  # Return both dendrogram and PCA results
  return(list(dendrogram = dend, pca = pca, pca_plot_file = pca_file))
}