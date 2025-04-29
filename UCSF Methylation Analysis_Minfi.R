#########################################
# UCSF Methylation Analysis Script
# Project: MPNST/pNF Methylation Profiling in NF1
# Author: ddzhang
# Date: 4/8/25
#########################################

# ----- 1. Load required libraries -----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# For methylation processing
if (!requireNamespace("minfi", quietly = TRUE))
  BiocManager::install("minfi")
if (!requireNamespace("IlluminaHumanMethylationEPICmanifest", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylationEPICmanifest")
if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE))
  BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# For heatmaps and clustering visualization
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
library(ggplot2)

if (!requireNamespace("cluster", quietly = TRUE))
  install.packages("cluster")
if (!requireNamespace("factoextra", quietly = TRUE))
  install.packages("factoextra")
library(cluster)
library(factoextra)

# ----- 2. Define Data Paths -----
# Directory where your raw IDAT files are stored:
idat_dir <- "/c4/home/ddzhang/MPNST_pNF_Methylation/UCSF/Raw_IDAT_Files/"

# Path to your clinical sample info (CSV) in your Methylation_data_reference folder:
metadata_path <- "/c4/home/ddzhang/Methylation_data_reference/SampleInfo_DNAM_20241114.csv"

# ----- 3. Load Clinical Metadata -----
sample_info <- read.csv(metadata_path, stringsAsFactors = FALSE)
head(sample_info)

# Create a shorter version of Basename (strip directory prefix)
# Adjust the column name "Basename" if your CSV uses a different name.
sample_info$Basename_short <- sub(".*/", "", sample_info$Basename)
# Set row names for easy matching later
rownames(sample_info) <- sample_info$Basename_short
# Inspect the modified metadata
head(sample_info)

# ----- 4. Load and Preprocess Methylation Data (Minfi Pipeline) -----
# Load the IDAT files from the UCSF directory.
RGSet <- read.metharray.exp(base = idat_dir, force = TRUE)
# Normalize using quantile normalization (preprocessQuantile, a lighter alternative to Noob)
MSet <- preprocessQuantile(RGSet)
# Extract beta values: rows = CpG sites, columns = samples
betas <- getBeta(MSet)
cat("Beta matrix dimensions (probes x samples):", dim(betas), "\n")

# Optional: Remove probes with any missing values
betas <- betas[complete.cases(betas), ]
cat("Beta matrix dimensions after filtering:", dim(betas), "\n")

# ----- 5. Exploratory Analysis: PCA & Hierarchical Clustering -----
# Transpose beta matrix so that rows represent samples.
betas_t <- t(betas)
# Run PCA on the transposed beta matrix.
pca_result <- prcomp(betas_t, scale. = TRUE)
# Basic PCA plot in base R
plot(pca_result$x[,1], pca_result$x[,2],
     main = "PCA of Methylation Data",
     xlab = "PC1", ylab = "PC2",
     pch = 19, col = "steelblue")

# Compute Euclidean distance among samples
dist_matrix <- dist(betas_t, method = "euclidean")
# Perform hierarchical clustering using Ward.D2 method
hc <- hclust(dist_matrix, method = "ward.D2")
# Plot dendrogram for clustering
plot(hc, main = "Hierarchical Clustering of Methylation Profiles", xlab = "", sub = "", cex = 0.8)

# ----- 6. Determine Optimal Number of Clusters Using Silhouette Analysis -----
# Test a few cluster numbers (k = 2 to 5)
for (k in 2:5) {
  clusters_temp <- cutree(hc, k = k)
  sil <- silhouette(clusters_temp, dist_matrix)
  avg_sil <- mean(sil[, 3])
  cat("Average silhouette width for k =", k, "is", avg_sil, "\n")
  # Optionally, plot silhouette
  plot(sil, main = paste("Silhouette Plot for k =", k))
}
# Based on your analysis and silhouette values, choose the optimal k:
optimal_k <- 2   # For example, if k=2 gives the best separation.
newClusters <- cutree(hc, k = optimal_k)
newClusters <- as.factor(newClusters)
cat("Final cluster assignment:\n")
print(table(newClusters))

# ----- 7. Merge New Cluster Assignment with Clinical Data -----
# Reorder sample_info to match the sample order in betas (using our short Basename)
sample_info_ordered <- sample_info[match(colnames(betas), sample_info$Basename_short), ]
# Add new cluster assignments as a new column (ensuring ordering is correct)
sample_info_ordered$NewCluster <- newClusters
head(sample_info_ordered)

# ----- 8. Generate an Annotated Heatmap -----
# We'll use a subset of probes for visualization, e.g., top 1000 variable probes:
top_var <- order(apply(betas, 1, var), decreasing = TRUE)[1:1000]

# For annotation, extract the "Tumor" and new cluster label from the metadata.
annotation_col <- sample_info_ordered[, c("Tumor", "NewCluster")]
# Optionally, rename columns for clarity:
colnames(annotation_col) <- c("Tumor", "Cluster")

# Plot the heatmap using pheatmap.
pheatmap(betas[top_var, ],
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Heatmap of Top 1000 Variable CpGs\nwith Tumor & Cluster Annotation")

# ----- 9. Create a Bar Plot of Cluster Distribution by Tumor Type -----
# Generate a contingency table comparing new clusters versus tumor type.
tbl <- table(sample_info_ordered$NewCluster, sample_info_ordered$Tumor)
df_tbl <- as.data.frame(tbl)
# Rename columns appropriately if needed.
colnames(df_tbl) <- c("Cluster", "Tumor", "Count")
# Plot with ggplot2
ggplot(df_tbl, aes(x = Cluster, y = Count, fill = Tumor)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Cluster",
       y = "Number of Samples",
       fill = "Tumor Type",
       title = "Distribution of Tumor Types by Cluster") +
  theme_minimal()

# ----- 10. (Optional) Save the Beta Matrix to CSV -----
write.csv(betas, file = "UCSF_minfi_beta_matrix.csv", row.names = TRUE)

# ----- End of Script -----
