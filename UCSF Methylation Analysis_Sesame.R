#########################################
# UCSF Sesame Methylation Analysis Script
# Project: MPNST/pNF Methylation Profiling in NF1
# Author: Daniel Zhang
# Date: 4/9/25
#########################################

# ----- 1. Load Required Libraries -----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("sesame", quietly = TRUE))
  BiocManager::install("sesame")
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("cluster", quietly = TRUE))
  install.packages("cluster")

library(sesame)
library(pheatmap)
library(ggplot2)
library(cluster)

# ----- 2. Define Data Paths -----
idat_dir <- "/c4/home/ddzhang/MPNST_pNF_Methylation/UCSF/Raw_IDAT_Files/"
metadata_path <- "/c4/home/ddzhang/Methylation_data_reference/SampleInfo_DNAM_20241114.csv"

# ----- 3. Load Clinical Metadata -----
sample_info <- read.csv(metadata_path, stringsAsFactors = FALSE)
# Create a shorter version of the Basename (adjust column name if necessary)
sample_info$Basename_short <- sub(".*/", "", sample_info$Basename)
rownames(sample_info) <- sample_info$Basename_short
print(head(sample_info))

# ----- 4. Identify Valid IDAT Basenames -----
# List all IDAT files in the directory
idat_files <- list.files(idat_dir, pattern = "idat$", full.names = TRUE)
# Get the unique basenames (removing _Grn.idat and _Red.idat suffixes)
idat_basenames <- unique(sub("_Grn\\.idat$|_Red\\.idat$", "", basename(idat_files)))
# Determine samples common to both the IDAT files and the metadata
valid_ids <- intersect(idat_basenames, sample_info$Basename_short)
cat("Number of valid samples:", length(valid_ids), "\n")

# ----- 5. Load and Preprocess IDAT Pairs with Progress Bar -----
pb <- txtProgressBar(min = 0, max = length(valid_ids), style = 3)
sdf_list <- lapply(seq_along(valid_ids), function(i) {
  setTxtProgressBar(pb, i)
  # Construct full file path by concatenating idat_dir and the valid ID
  readIDATpair(paste0(idat_dir, valid_ids[i]))
})
close(pb)

# ----- 6. Extract Beta Values -----
system.time({
  beta_list <- lapply(sdf_list, getBetas)
})
# Combine beta_list into a matrix (rows: CpG probes, columns: samples)
betas <- do.call(cbind, beta_list)
colnames(betas) <- valid_ids  # Use valid_ids as column names
# Optionally, remove probes with any missing values:
betas <- betas[complete.cases(betas), ]
cat("Beta matrix dimensions (probes x samples):", dim(betas), "\n")

# ----- 6B. Remove sex chromosome probes (X and Y) -----
probe_annot <- sesameDataGet("probe.features")  # load probe annotations
autosomal_probes <- rownames(probe_annot)[probe_annot$chr %in% paste0("chr", 1:22)]
# Filter beta matrix to keep only autosomal probes
betas <- betas[rownames(betas) %in% autosomal_probes, ]
cat("Beta matrix after removing sex chromosome probes:", dim(betas), "\n")


# ----- 7. Reorder and Subset Metadata to Match Beta Matrix -----
# Subset sample_info to include only samples in valid_ids
sample_info <- sample_info[sample_info$Basename_short %in% valid_ids, ]
# Reorder sample_info so that its row order matches betas column order
sample_info <- sample_info[match(colnames(betas), sample_info$Basename_short), ]
# Confirm alignment
if (all(colnames(betas) == sample_info$Basename_short)) {
  cat("Sample metadata and beta matrix column names are aligned.\n")
} else {
  stop("Mismatch between beta matrix column names and sample_info!")
}

# ----- 8. Exploratory PCA -----
# Transpose beta matrix so that each row represents a sample.
betas_t <- t(betas)
# Run principal component analysis with scaling.
pca_result <- prcomp(betas_t, scale. = TRUE)
# Basic PCA plot for PC1 vs. PC2
plot(
  pca_result$x[,1], pca_result$x[,2],
  main = "PCA of Sesame Methylation Data",
  xlab = "PC1", ylab = "PC2",
  pch = 19, col = "steelblue"
)

# ----- 9. Hierarchical Clustering and Dendrogram -----
dist_matrix <- dist(betas_t, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")
plot(hc, main = "Hierarchical Clustering of Sesame Profiles", xlab = "", sub = "", cex = 0.8)

# ----- 10. Silhouette Analysis for Optimal Cluster Determination -----
for (k in 2:5) {
  clusters_temp <- cutree(hc, k = k)
  sil <- silhouette(clusters_temp, dist_matrix)
  avg_sil <- mean(sil[,3])
  cat("Average silhouette width for k =", k, ":", avg_sil, "\n")
  # Optionally, visualize silhouette plot:
  # plot(sil, main = paste("Silhouette Plot for k =", k))
}
# Choose the optimal number of clusters (adjust based on silhouette output)
optimal_k <- 2  
final_clusters <- factor(cutree(hc, k = optimal_k))

# ----- 11. Merge Cluster Assignment with Metadata -----
sample_info$NewCluster <- final_clusters
print(head(sample_info))

# ----- 12. Identify Top Variable Probes for Visualization -----
# Calculate variance for each probe and select top 1000 most variable ones.
var_probes <- apply(betas, 1, var)
top_var <- order(var_probes, decreasing = TRUE)[1:1000]

# ----- 13. Generate Heatmap Annotated by Cluster -----
annotation_cluster <- data.frame(Cluster = sample_info$NewCluster)
rownames(annotation_cluster) <- colnames(betas)
pheatmap(
  betas[top_var, ],
  show_rownames = FALSE,
  annotation_col = annotation_cluster,
  main = "Heatmap of Top 1000 Variable Probes\n(Annotated by Cluster)"
)

# ----- 14. Generate Heatmap Annotated by Tumor Type -----
# (Assumes that your sample_info contains a 'Tumor' column)
annotation_tumor <- data.frame(Tumor = sample_info$Tumor)
rownames(annotation_tumor) <- colnames(betas)
pheatmap(
  betas[top_var, ],
  show_rownames = FALSE,
  annotation_col = annotation_tumor,
  main = "Heatmap of Top 1000 Variable Probes\n(Annotated by Tumor Type)"
)

# ----- 15. Generate Bar Plot of Cluster Distribution by Tumor Type -----
tbl <- table(sample_info$NewCluster, sample_info$Tumor)
df_tbl <- as.data.frame(tbl)
colnames(df_tbl) <- c("Cluster", "Tumor", "Count")
ggplot(df_tbl, aes(x = Cluster, y = Count, fill = Tumor)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = "Cluster",
    y = "Number of Samples",
    fill = "Tumor Type",
    title = "Distribution of Tumor Types by Cluster (Sesame)"
  ) +
  theme_minimal()

# ----- 16. (Optional) Export the Beta Matrix to CSV -----
write.csv(betas, file = "UCSF_sesame_beta_matrix.csv", row.names = TRUE)

#########################################
# End of UCSF Sesame Methylation Analysis Script
#########################################
