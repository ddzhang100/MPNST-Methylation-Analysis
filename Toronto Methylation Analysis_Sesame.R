#########################################
# Toronto Sesame Methylation Analysis Script
# Using ConsensusClusterPlus with Dendrogram-Aligned Heatmap
# Project: MPNST/pNF Methylation Profiling in NF1 (Toronto Cohort)
# Author: [Your Name]
# Date: [Current Date]
#########################################

# ----- 1. Load Required Libraries -----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Sesame and associated packages
if (!requireNamespace("sesame", quietly = TRUE))
  BiocManager::install("sesame")
if (!requireNamespace("sesameData", quietly = TRUE))
  BiocManager::install("sesameData")
library(sesame)

# Visualization, clustering and consensus clustering
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("cluster", quietly = TRUE))
  install.packages("cluster")
if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE))
  BiocManager::install("ConsensusClusterPlus")
library(pheatmap)
library(ggplot2)
library(cluster)
library(ConsensusClusterPlus)

# For reading Excel files
if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)

# ----- 2. Define File Paths -----
# Path to unzipped IDAT files (Toronto cohort)
idat_dir <- "/c4/home/ddzhang/MPNST_pNF_Methylation/Toronto/GSE207207_RAW"
# Path to the metadata Excel file (reduced metadata from GEO series matrix)
meta_reduced_path <- "/c4/home/ddzhang/Methylation_data_reference/GSE207207_series_matrix_reduced.xlsx"

# ----- 3. Load Toronto Metadata from Excel -----
meta_reduced <- read.xlsx(meta_reduced_path, sheet = 1)
# Assume the metadata contains at least columns: geo_accession and tumor_type
rownames(meta_reduced) <- meta_reduced$geo_accession
cat("Distinct tumor types in metadata:\n")
print(unique(meta_reduced$tumor_type))

# ----- 4. Process IDAT Filenames and Map to Metadata -----
# List all .idat files in the directory (including subfolders)
idat_files <- list.files(path = idat_dir, pattern = "idat$", full.names = TRUE, recursive = TRUE)
cat("Total number of raw IDAT files found:", length(idat_files), "\n")

# Remove the _Grn.idat or _Red.idat suffix to obtain the full sample prefix
prefixes <- unique(sub("(_Grn|_Red)\\.idat$", "", basename(idat_files)))
cat("Unique full sample prefixes found (first 5):\n")
print(head(prefixes, 5))

# Extract the short GSM ID from each prefix (assume it is everything before the first underscore)
sample_ids <- sub("_.*", "", prefixes)
# Create a mapping data frame between full prefix and the GSM ID
prefix_map <- data.frame(prefix = prefixes, gsm = sample_ids, stringsAsFactors = FALSE)
cat("Prefix mapping (first 5 rows):\n")
print(head(prefix_map, 5))

# Filter the mapping to only samples that appear in meta_reduced based on the GSM ID
common_gsm <- intersect(prefix_map$gsm, rownames(meta_reduced))
prefix_map_filtered <- prefix_map[prefix_map$gsm %in% common_gsm, ]
cat("Number of samples with matching metadata:", nrow(prefix_map_filtered), "\n")

# ----- 5. Load IDAT Files Using Sesame and Normalize -----
# Time the IDAT loading process
start_time <- Sys.time()
sdf_list <- lapply(prefix_map_filtered$prefix, function(px) {
  # Full path to sample prefix; assumes that files are in the top-level idat_dir
  full_path <- file.path(idat_dir, px)
  readIDATpair(full_path)  # This looks for full_path_Grn.idat and full_path_Red.idat
})
end_time <- Sys.time()
cat("Time taken to load IDAT pairs:", end_time - start_time, "\n")

# Normalize the data for each sample using noob and dye bias correction
sdf_norm <- lapply(sdf_list, function(sdf) {
  sdf <- noob(sdf)
  sdf <- dyeBiasCorrTypeINorm(sdf)
  sdf
})
# Extract beta values and combine into a matrix (rows = probes, columns = samples)
beta_list <- lapply(sdf_norm, getBetas)
betas <- do.call(cbind, beta_list)
# Set column names to the short GSM IDs
colnames(betas) <- prefix_map_filtered$gsm
# Remove any probes with missing values
betas <- betas[complete.cases(betas), ]
cat("Beta matrix dimensions after filtering (probes x samples):", dim(betas), "\n")

# ----- 6. Align Beta Matrix with Metadata -----
common_ids_final <- intersect(colnames(betas), rownames(meta_reduced))
betas <- betas[, common_ids_final, drop = FALSE]
meta_reduced <- meta_reduced[common_ids_final, , drop = FALSE]
stopifnot(all(colnames(betas) == rownames(meta_reduced)))
cat("Final number of aligned samples:", ncol(betas), "\n")

# ----- 7. Exploratory PCA -----
betas_t <- t(betas)
pca_result <- prcomp(betas_t, scale. = TRUE)
# Basic PC1 vs. PC2 plot
plot(pca_result$x[,1], pca_result$x[,2],
     main = "PCA of Sesame Beta Matrix (Toronto)",
     xlab = paste("PC1 (", round(100 * (pca_result$sdev[1]^2) / sum(pca_result$sdev^2), 2), "%)", sep = ""),
     ylab = paste("PC2 (", round(100 * (pca_result$sdev[2]^2) / sum(pca_result$sdev^2), 2), "%)", sep = ""),
     pch = 19, col = "steelblue")

# Additional PCA plots: PC1 vs. PC3 and PC2 vs. PC3
pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 2)
# PC1 vs PC3:
plot(pca_result$x[,1], pca_result$x[,3],
     main = "PCA: PC1 vs PC3",
     xlab = paste("PC1 (", pca_var_perc[1], "%)", sep = ""),
     ylab = paste("PC3 (", pca_var_perc[3], "%)", sep = ""),
     pch = 19, col = "steelblue")
# PC2 vs PC3:
plot(pca_result$x[,2], pca_result$x[,3],
     main = "PCA: PC2 vs PC3",
     xlab = paste("PC2 (", pca_var_perc[2], "%)", sep = ""),
     ylab = paste("PC3 (", pca_var_perc[3], "%)", sep = ""),
     pch = 19, col = "steelblue")

# ----- 8. Hierarchical Clustering -----
dist_matrix <- dist(betas_t, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")
plot(hc, main = "Hierarchical Clustering (Toronto - Sesame)",
     xlab = "", sub = "", cex = 0.8)

# ----- 9. Consensus Clustering for Robust Cluster Assignment -----
# Select the top 1000 variable probes (for consensus clustering)
var_probes_all <- apply(betas, 1, var)
top_var_indices <- order(var_probes_all, decreasing = TRUE)[1:1000]
betas_top_subset <- betas[top_var_indices, , drop = FALSE]

# Run ConsensusClusterPlus on the top variable probes
cc_results <- ConsensusClusterPlus(as.matrix(betas_top_subset),
                                   maxK = 6,
                                   reps = 500,
                                   pItem = 0.8,
                                   pFeature = 1,
                                   clusterAlg = "hc",
                                   distance = "euclidean",
                                   seed = 1234,
                                   title = "Consensus_Results_Toronto_Sesame",
                                   plot = "png")
# Suppose k = 2 is optimal based on CDF and delta area plots
consensus_clusters <- cc_results[[2]]$consensusClass
# Ensure cluster assignments match the order of betas columns
consensus_clusters <- consensus_clusters[match(colnames(betas), names(consensus_clusters))]
meta_reduced$NewCluster <- factor(consensus_clusters)
cat("Consensus clustering assignments added to metadata.\n")

# ----- 10. Forcing Dendrogram Matching for the Heatmap -----
# Compute column dendrogram for the top variable probes
dist_cols <- dist(t(betas_top_subset), method = "euclidean")
hc_cols <- hclust(dist_cols, method = "ward.D2")
plot(hc_cols, main = "Column Dendrogram (Top 1000 Probes) - Toronto (Sesame)",
     xlab = "", sub = "", cex = 0.8)

# Compute row dendrogram on the same set
dist_rows <- dist(betas_top_subset, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "ward.D2")

# ----- 11. Generate Heatmap Annotated by Tumor Type and Consensus Cluster -----
# Build annotation data frame from metadata (ensure row names match betas columns)
annotation_col <- data.frame(Tumor = meta_reduced$tumor_type,
                             Cluster = meta_reduced$NewCluster)
rownames(annotation_col) <- rownames(meta_reduced)

pheatmap(
  mat = betas_top_subset,
  cluster_cols = hc_cols,     # Force using our precomputed column dendrogram
  cluster_rows = hc_rows,
  annotation_col = annotation_col,
  treeheight_col = 100,       # Increase dendrogram height for clarity
  show_rownames = FALSE,
  main = "Heatmap of Top 1000 Variable Probes\n(Precomputed Column Dendrogram + Tumor Annotation - Toronto - Sesame)"
)

# ----- 12. Generate Bar Plot of Consensus Cluster Distribution by Tumor Type -----
tbl <- table(meta_reduced$NewCluster, meta_reduced$tumor_type)
df_tbl <- as.data.frame(tbl)
colnames(df_tbl) <- c("Cluster", "TumorType", "Count")
ggplot(df_tbl, aes(x = Cluster, y = Count, fill = TumorType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Consensus Cluster", y = "Number of Samples",
       fill = "Tumor Type",
       title = "Distribution of Tumor Types by Consensus Cluster (Toronto - Sesame)") +
  theme_minimal()

# ----- 13. (Optional) Export Beta Matrix to CSV -----
write.csv(betas, file = "Toronto_sesame_beta_matrix.csv", row.names = TRUE)

#########################################
# End of Toronto Sesame Methylation Analysis Script
#########################################
