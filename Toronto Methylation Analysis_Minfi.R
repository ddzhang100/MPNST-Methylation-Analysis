#########################################
# Toronto Minfi Methylation Analysis with Tumor Type Annotation
# Using meta_reduced for tumor info, consensus clustering, and dendrogram + heatmap
# Project: MPNST/pNF Methylation Profiling in NF1
# Author: [Your Name]
# Date: [Current Date]
#########################################

# ----- 1. Load Required Libraries -----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Minfi and annotation packages
if (!requireNamespace("minfi", quietly = TRUE))
  BiocManager::install("minfi")
if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

# Visualization, clustering, and consensus clustering packages
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("cluster", quietly = TRUE))
  install.packages("cluster")
if (!requireNamespace("factoextra", quietly = TRUE))
  install.packages("factoextra")
if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE))
  BiocManager::install("ConsensusClusterPlus")

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(pheatmap)
library(ggplot2)
library(cluster)
library(factoextra)
library(ConsensusClusterPlus)

# ----- 2. Define Data Paths -----
# Path to your unzipped IDAT files for the Toronto cohort
idat_dir <- "/c4/home/ddzhang/MPNST_pNF_Methylation/Toronto/GSE207207_RAW"

# The meta_reduced CSV that you created with 'geo_accession' and 'tumor_type'
meta_reduced_path <- "/c4/home/ddzhang/Methylation_data_reference/GSE207207_series_matrix_reduced.csv"

# ----- 3. Load the Methylation Data (Minfi Pipeline) -----
RGSet <- read.metharray.exp(base = idat_dir, force = TRUE, verbose = TRUE)
MSet <- preprocessQuantile(RGSet)
betas <- getBeta(MSet)
cat("Beta matrix dimensions (probes x samples):", dim(betas), "\n")

# Remove probes with any missing data
betas <- betas[complete.cases(betas), ]
cat("Beta matrix dimensions after filtering:", dim(betas), "\n")

# ----- 4. Load Tumor Metadata from meta_reduced -----
# This file should have columns: geo_accession, tumor_type
meta_reduced <- read.csv(meta_reduced_path, stringsAsFactors = FALSE)
# Set rownames to the geo_accession for easy matching
rownames(meta_reduced) <- meta_reduced$geo_accession

# Inspect first few rows of metadata
head(meta_reduced)

# ----- 5. Align Beta Matrix Columns with Metadata Rows -----
# If the column names of betas are something like "GSM6281328_2022A209716_R07C01",
# but meta_reduced$geo_accession is just "GSM6281328", we need to strip the IDAT suffix.
# For example:
colnames(betas) <- sub("_.*", "", colnames(betas))  # remove everything after first underscore
# Now colnames(betas) should match "GSM6281328", etc.

# Filter and reorder betas + metadata to common IDs
common_ids <- intersect(colnames(betas), rownames(meta_reduced))
betas <- betas[, common_ids, drop = FALSE]
meta_reduced <- meta_reduced[common_ids, , drop = FALSE]
stopifnot(all(colnames(betas) == rownames(meta_reduced)))

cat("Number of samples retained after alignment:", ncol(betas), "\n")

# ----- 6. Exploratory PCA -----
betas_t <- t(betas)
pca_result <- prcomp(betas_t, scale. = TRUE)
plot(pca_result$x[,1], pca_result$x[,2],
     main = "PCA of Methylation Data (Toronto - Minfi)",
     xlab = "PC1", ylab = "PC2",
     pch = 19, col = "steelblue")

# ----- 7. Hierarchical Clustering and Dendrogram -----
dist_matrix <- dist(betas_t, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")
plot(hc, main = "Hierarchical Clustering (Toronto - Minfi)",
     xlab = "", sub = "", cex = 0.8)

# ----- 8. Consensus Clustering for Robust Cluster Assignment -----
# Identify top 1000 variable probes
var_probes <- apply(betas, 1, var)
top_var_indices <- order(var_probes, decreasing = TRUE)[1:1000]
betas_top <- betas[top_var_indices, , drop = FALSE]

# Run consensus clustering
cc_results <- ConsensusClusterPlus(as.matrix(betas_top),
                                   maxK = 6,
                                   reps = 500,
                                   pItem = 0.8,
                                   pFeature = 1,
                                   clusterAlg = "hc",
                                   distance = "euclidean",
                                   seed = 1234,
                                   title = "Consensus_Results_Toronto",
                                   plot = "png")
# Suppose we pick k=2 based on the CDF/delta area plot
consensus_clusters <- cc_results[[2]]$consensusClass
# Match cluster assignments to columns of betas
consensus_clusters <- consensus_clusters[match(colnames(betas), names(consensus_clusters))]
# Store in metadata
meta_reduced$NewCluster <- factor(consensus_clusters)

cat("Consensus clusters added to meta_reduced.\n")

# ----- 9. Forced Column Ordering for Heatmap & Precomputed Dendrogram -----
dist_cols <- dist(t(betas_top), method = "euclidean")
hc_cols <- hclust(dist_cols, method = "ward.D2")

# Plot the new column dendrogram for reference
plot(hc_cols, main = "Column Dendrogram (Top 1000 Probes) - Toronto",
     xlab = "", sub = "", cex = 0.8)

# row clustering
dist_rows <- dist(betas_top, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "ward.D2")

# ----- 10. Generate Heatmap Annotated by Tumor Type -----
# Create an annotation data frame for the heatmap
annotation_col <- data.frame(Tumor = meta_reduced$tumor_type,
                             Cluster = meta_reduced$NewCluster)
rownames(annotation_col) <- rownames(meta_reduced)

# Option A: Show the top dendrogram by passing hc_cols
pheatmap(
  mat = betas_top,
  cluster_cols = hc_cols,
  cluster_rows = hc_rows,
  annotation_col = annotation_col,
  treeheight_col = 100,
  show_rownames = FALSE,
  main = "Heatmap of Top 1000 Variable Probes - Toronto\nPrecomputed Dendrogram + Tumor Annotation"
)

# ----- 11. (Optional) Bar Plot of Tumor Type vs. Clusters -----
tbl <- table(meta_reduced$NewCluster, meta_reduced$tumor_type)
df_tbl <- as.data.frame(tbl)
colnames(df_tbl) <- c("Cluster", "TumorType", "Count")
ggplot(df_tbl, aes(x = Cluster, y = Count, fill = TumorType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Consensus Cluster", y = "Number of Samples",
       fill = "Tumor Type",
       title = "Distribution of Tumor Types by Consensus Cluster (Toronto)") +
  theme_minimal()

# ----- 12. (Optional) Export the Beta Matrix to CSV -----
write.csv(betas, file = "Toronto_minfi_beta_matrix.csv", row.names = TRUE)

#########################################
# End of Toronto Minfi Methylation Analysis Script with Tumor Type
#########################################
