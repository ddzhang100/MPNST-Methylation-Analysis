##############################################
# MPNST vs. Normal‑nerve‑sheath (NF) DMRcate Pipeline
# – From raw IDATs to region‑level DMRs and enrichment
# Author: Daniel Zhang
# Date: 2025‑04‑16
##############################################

### 1. Load & install required packages ###
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

pkgs <- c(
  "minfi",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "DMRcate",
  "limma",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ggplot2"
)
for (p in pkgs) {
  if (!requireNamespace(p, quietly=TRUE))
    BiocManager::install(p)
  library(p, character.only=TRUE)
}

### 2. Define paths & read metadata ###
idat_dir      <- "/c4/home/ddzhang/MPNST_pNF_Methylation/UCSF/Raw_IDAT_Files/"
metadata_path <- "/c4/home/ddzhang/Methylation_data_reference/SampleInfo_DNAM_20241114.csv"

sample_info <- read.csv(metadata_path, stringsAsFactors=FALSE)
sample_info$Basename_short <- sub(".*/", "", sample_info$Basename)
rownames(sample_info) <- sample_info$Basename_short
stopifnot("Tumor" %in% colnames(sample_info))  # must be "MPNST" or "NF"
cat("Metadata loaded for", nrow(sample_info), "samples.\n")

### 3. Read & normalize raw data (minfi) ###
RGSet <- read.metharray.exp(base=idat_dir, force=TRUE)
cat("Loaded", length(RGSet), "IDAT samples.\n")

MSet  <- preprocessQuantile(RGSet)
betas <- getBeta(MSet)
cat("Raw beta matrix dims:", dim(betas), "\n")

# filter out any probes with missing values
betas <- betas[complete.cases(betas), ]
cat("Filtered beta dims:", dim(betas), "\n")

# Save for reproducibility
save(betas, file="betas_minfi_common.RData")

### 4. Convert to M‑values ###
M_vals <- log2(betas / (1 - betas))
M_vals[!is.finite(M_vals)] <- NA
cat("M‑value matrix dims:", dim(M_vals), "\n")

### 5. Align samples between M‑values and metadata ###
common_samples <- intersect(colnames(M_vals), rownames(sample_info))
M_vals         <- M_vals[, common_samples, drop=FALSE]
sample_info    <- sample_info[common_samples, , drop=FALSE]
cat("Using", length(common_samples), "common samples.\n")

### 6. Build design & contrast for MPNST vs NF ###
group_factor <- factor(sample_info$Tumor, levels=c("NF","MPNST"))
design_tumor <- model.matrix(~0 + group_factor)
colnames(design_tumor) <- c("NF","MPNST")
cat("Design dims:", dim(design_tumor), "\n")
print(colnames(design_tumor))

contrast_tumor <- makeContrasts(
  Tumor_vs_NF = MPNST - NF,
  levels      = design_tumor
)
print(contrast_tumor)

### 7. Annotate CpGs & call DMRs (DMRcate) ###
annot_tumor <- cpg.annotate(
  object        = M_vals,
  datatype      = "array",
  what          = "M",
  analysis.type = "differential",
  design        = design_tumor,
  contrasts     = TRUE,
  cont.matrix   = contrast_tumor,
  coef          = "Tumor_vs_NF",
  fdr           = 0.05,
  arraytype     = "EPIC"    # EPIC v1; use "450K" if 450K array
)
cat("Probe‑wise annotation complete.\n")

dmrc_tumor <- dmrcate(annot_tumor, lambda=1000, C=2)
cat("DMRcate smoothing complete.\n")

dmr_ranges_tumor <- extractRanges(dmrc_tumor, genome="hg19")
cat("Number of DMRs (MPNST vs NF):", length(dmr_ranges_tumor), "\n")

# save DMR table
dmr_df_tumor <- as.data.frame(dmr_ranges_tumor)
write.csv(dmr_df_tumor,
          file="DMRs_MPNST_vs_NF.csv",
          row.names=FALSE)
cat("DMR results written to DMRs_MPNST_vs_NF.csv\n")

### 8. Extract DMR‑associated genes ###
gene_strings_tumor <- mcols(dmr_ranges_tumor)$overlapping.genes
gene_list_tumor   <- unique(
  unlist(strsplit(as.character(gene_strings_tumor), split="[,;]"))
)
gene_list_tumor   <- trimws(gene_list_tumor[nzchar(gene_list_tumor)])
cat("Unique DMR genes:", length(gene_list_tumor), "\n")

### 9. GO (BP) enrichment ###
ego_tumor <- enrichGO(
  gene          = gene_list_tumor,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
cat("Top GO (BP) terms:\n")
print(head(as.data.frame(ego_tumor), 10))

# plot
barplot(ego_tumor, showCategory=20,
        title="GO (BP) Enrichment: MPNST vs NF")

### 10. KEGG enrichment ###
library(clusterProfiler)  # ensure loaded
ids_tumor <- bitr(
  gene_list_tumor,
  fromType="SYMBOL",
  toType="ENTREZID",
  OrgDb=org.Hs.eg.db
)
ekegg_tumor <- enrichKEGG(
  gene         = ids_tumor$ENTREZID,
  organism     = 'hsa',
  pAdjustMethod= "BH",
  qvalueCutoff = 0.05
)
cat("Top KEGG pathways:\n")
print(head(as.data.frame(ekegg_tumor), 10))

# plot
barplot(ekegg_tumor, showCategory=20,
        title="KEGG Enrichment: MPNST vs NF")

##############################################
# End of MPNST vs NF DMRcate + Enrichment Pipeline
##############################################
