#!/usr/bin/env Rscript

# ============================================================
# Minimal, reproducible pipeline retained from the original scripts:
#   Raw 10x Xenium/CellRanger count matrix
#     -> Seurat QC + clustering (for CM cluster selection)
#     -> sgRNA assignment (single-KO cells)
#     -> perturbation-signature object (CalcPerturbSig; optional but kept)
#     -> module scores per KO
#     -> sgNT (NTgRNA1) vs each sgRNA comparison (Wilcoxon test)
#
# Outputs (written under ./results/):
#   - rds/PS_peturb.rds
#   - tables/module_score_long.csv
#   - tables/module_score_stats_vs_NT.csv
#   - plots/module_scores_violin_by_KO.pdf
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(scales)
  library(readxl)
})

# -----------------------------
# 0) Project paths (EDIT THESE)
# -----------------------------

# Use the folder containing this script as the "project root" by default.
# If you run interactively, setwd() to the project root first.
PROJECT_ROOT <- getwd()

# 10x filtered_feature_bc_matrix folder
TENX_DIR <- file.path(PROJECT_ROOT, "cellranger/outs/count/filtered_feature_bc_matrix")

# Optional: a probe-set CSV used in the original analysis (kept for traceability)
# If you don't have it, you can comment it out.
PROBESET_CSV <- file.path(PROJECT_ROOT, "cellranger/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10_custom.csv")

# Mapping file: sgRNA ID -> gene name (as in original script)
SGRNA2GENE_CSV <- file.path(PROJECT_ROOT, "rObjs/sgRNA2gene.csv")

# Excel sheet defining maturation-related gene modules (Category / gene columns)
MODULE_XLSX <- file.path(PROJECT_ROOT, "Target.gene.100523.xlsx")

# Output dirs
OUT_RDS   <- file.path(PROJECT_ROOT, "results/rds")
OUT_PLOTS <- file.path(PROJECT_ROOT, "results/plots")
OUT_TABS  <- file.path(PROJECT_ROOT, "results/tables")
dir.create(OUT_RDS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TABS, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) Load helpers
# -----------------------------
source(file.path(PROJECT_ROOT, "ps_functions_min.R"))

# -----------------------------
# 2) Read raw 10x -> Seurat
# -----------------------------
PS.data <- Read10X(data.dir = TENX_DIR)
PS.lib  <- CreateSeuratObject(counts = PS.data, project = "PS_lib")

# Mito percentage (match original: "^mt-" for mouse panel; switch to "^MT-" for human)
PS.lib[["percent.mt"]] <- PercentageFeatureSet(PS.lib, pattern = "^mt-")

# Robust outlier filter using MAD (kept from original)
is_outlier <- function(seu, metric, nmads) {
  M <- seu@meta.data[[metric]]
  med <- median(M, na.rm = TRUE)
  mad_value <- mad(M, constant = 1, na.rm = TRUE)
  lower <- med - nmads * mad_value
  upper <- med + nmads * mad_value
  (M < lower) | (M > upper)
}

PS.lib[["count_outlier"]] <- is_outlier(PS.lib, "nCount_RNA", 5) | is_outlier(PS.lib, "nFeature_RNA", 5)
PS.lib[["mt_outlier"]]    <- is_outlier(PS.lib, "percent.mt", 5)
PS.lib.filtered <- subset(PS.lib, subset = count_outlier == FALSE & mt_outlier == FALSE)

# Basic QC plots
pdf(file.path(OUT_PLOTS, "QC_pre_vs_post_filter.pdf"), onefile = TRUE, width = 12, height = 5)
print(VlnPlot(PS.lib.filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
print(VlnPlot(PS.lib,          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
print(FeatureScatter(PS.lib, feature1 = "nCount_RNA", feature2 = "percent.mt") +
        FeatureScatter(PS.lib.filtered, feature1 = "nCount_RNA", feature2 = "percent.mt"))
dev.off()

# -----------------------------
# 3) Normalize + cluster (for CM selection)
# -----------------------------
PS.lib.filtered <- NormalizeData(PS.lib.filtered)
PS.lib.filtered <- FindVariableFeatures(PS.lib.filtered, selection.method = "vst", nfeatures = 2000)
PS.lib.filtered <- ScaleData(PS.lib.filtered, features = rownames(PS.lib.filtered))
PS.lib.filtered <- RunPCA(PS.lib.filtered, features = VariableFeatures(PS.lib.filtered))

# Choose PCs using the same heuristic as original
pct  <- PS.lib.filtered[["pca"]]@stdev / sum(PS.lib.filtered[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1  <- which(cumu > 90 & pct < 5)[1]
co2  <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
pcs  <- min(co1, co2)

PS.lib.filtered <- FindNeighbors(PS.lib.filtered, dims = 1:pcs)
PS.lib.filtered <- FindClusters(PS.lib.filtered, resolution = 0.8)
PS.lib.filtered <- RunUMAP(PS.lib.filtered, dims = 1:pcs)

pdf(file.path(OUT_PLOTS, "UMAP_clusters.pdf"), width = 10, height = 6)
print(DimPlot(PS.lib.filtered, reduction = "umap", label = TRUE))
dev.off()

# -----------------------------
# 4) Subset cardiomyocytes (CM) by clusters
# -----------------------------
# IMPORTANT: This is the *exact* cluster selection used in the original script.
# Please verify cluster identities using marker plots if you apply this to a new dataset.
#
# Original note: "cluster 14,12,17, and 8 are removed as non-myocyte cluster"
# Code used: subset(..., idents = c(9,12,14,15,16,17), invert = T)
#
PS.lib.CM <- subset(PS.lib.filtered, idents = c(9, 12, 14, 15, 16, 17), invert = TRUE)

# -----------------------------
# 5) Extract sgRNA features + assign "true_sgRNA" for single-KO cells
# -----------------------------
# The original analysis extracted sgRNA rows by fixed indices in the count matrix.
# Those indices are specific to the panel/design; keep them only if your matrix layout matches.
# If your features are named consistently, a safer approach is to grep by prefix and select rows by name.
#
# Here we keep the original behavior for reproducibility.
PS.rna.df <- as.data.frame(PS.lib.CM@assays$RNA@layers$counts)
colnames(PS.rna.df) <- colnames(PS.lib.CM)
rownames(PS.rna.df) <- rownames(PS.lib.CM)

# (Optional) read probe set file (not required for the downstream module-score steps)
if (file.exists(PROBESET_CSV)) {
  custome.probe <- read.csv(PROBESET_CSV, skip = 55952, header = FALSE)
}

# Original hard-coded slicing (VERIFY for your panel):
PS.grna.df <- PS.rna.df[19048:19153, , drop = FALSE]

# Save a copy for transparency
write.csv(PS.grna.df, file.path(PROJECT_ROOT, "results/tables/gRNA_df_rawcounts.csv"), quote = FALSE)

# Build "double_expressing_cells" list per guide-pair (two rows per target)
double_expressing_cells <- list()
for (i in 1:(nrow(PS.grna.df) / 2)) {
  idx1 <- (i * 2) - 1
  idx2 <- i * 2
  gRNA1 <- rownames(PS.grna.df)[idx1]
  gRNA2 <- rownames(PS.grna.df)[idx2]
  gRNA_name <- rownames(PS.grna.df)[idx1]
  double_expressing_cells[[gRNA_name]] <- get_double_expressing_cells(PS.lib.CM, gRNA1, gRNA2, threshold = 2)
}

# Construct metadata table for sgRNA assignment
df <- data.frame(row_name = colnames(PS.lib.CM))
for (nm in names(double_expressing_cells)) {
  df[[nm]] <- df$row_name %in% double_expressing_cells[[nm]]
}
rownames(df) <- df$row_name
df <- df[, -1, drop = FALSE]

# Keep only cells with exactly one TRUE entry across targets (single-KO assignment)
subset_df <- df[rowSums(df) == 1, , drop = FALSE]
subset_df$true_sgRNA <- apply(subset_df, 1, get_true_sgRNA_names)

# Subset Seurat object and attach metadata
PS.lib.CM.assigned.singleKO <- subset(PS.lib.CM, cells = rownames(subset_df))
PS.lib.CM.assigned.singleKO <- AddMetaData(PS.lib.CM.assigned.singleKO, subset_df)

# Quick QC plot: number of cells per KO label
pdf(file.path(OUT_PLOTS, "KO_cell_barplot.pdf"), width = 12, height = 4)
tmp_meta <- data.frame(true_sgRNA = PS.lib.CM.assigned.singleKO$true_sgRNA)
print(
  ggplot(tmp_meta, aes(x = true_sgRNA)) +
    geom_bar() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "KO label (true_sgRNA)", y = "Cell count")
)
dev.off()

# -----------------------------
# 6) (Optional but kept) Build perturbation signature assay (PRTB)
# -----------------------------
# NOTE: This requires Seurat's perturb-seq utilities (CalcPerturbSig).
# If CalcPerturbSig is not available in your Seurat version, you can skip this block
# and compute module scores directly on PS.lib.CM.assigned.singleKO instead.
PS.lib.CM.assigned.singleKO <- FindVariableFeatures(PS.lib.CM.assigned.singleKO, selection.method = "vst", nfeatures = 2000)
PS.lib.CM.assigned.singleKO <- ScaleData(PS.lib.CM.assigned.singleKO, features = rownames(PS.lib.CM.assigned.singleKO))
PS.lib.CM.assigned.singleKO <- RunPCA(PS.lib.CM.assigned.singleKO, features = VariableFeatures(PS.lib.CM.assigned.singleKO))

PS_peturb <- CalcPerturbSig(
  object = PS.lib.CM.assigned.singleKO,
  assay = "RNA",
  slot = "data",
  gd.class = "true_sgRNA",
  nt.cell.class = "NTgRNA1",
  reduction = "pca",
  ndims = 40,
  num.neighbors = 20,
  new.assay.name = "PRTB"
)

# Save the object for downstream reuse
saveRDS(PS_peturb, file.path(OUT_RDS, "PS_peturb.rds"))

# -----------------------------
# 7) Compute module scores on PS_peturb
# -----------------------------
# Read module gene list Excel (expects columns: Category, gene)
m.gene <- readxl::read_excel(MODULE_XLSX)

stopifnot(all(c("Category", "gene") %in% colnames(m.gene)))

cd_features <- list()
for (cat in unique(m.gene$Category)) {
  cd_features[[cat]] <- dplyr::filter(m.gene, Category == cat)$gene
}

# Select the modules you want to score (kept from original script)
select_module <- c(
  "Transcription factor", "Fatty acid", "T-tubule", "Mitochondria",
  "Channel and handling", "Alternative splicing gene", "Myofilament",
  "Growth", "Kinase", "RNA splicing factor"
)

# Add module scores into metadata
for (key in select_module) {
  if (!key %in% names(cd_features)) {
    warning(sprintf("Module '%s' not found in MODULE_XLSX; skipping.", key))
    next
  }
  mod_name <- paste0(gsub(" ", "_", key), "_module_score")
  PS_peturb[[mod_name]] <- module_score(cd_features[[key]], seurat = PS_peturb, assay = "RNA")
}

# -----------------------------
# 8) sgNT vs each sgRNA comparison (Wilcoxon)
# -----------------------------
# Create long table of scores
score_cols <- grep("_module_score$", colnames(PS_peturb@meta.data), value = TRUE)
stopifnot("true_sgRNA" %in% colnames(PS_peturb@meta.data))

scores_long <- PS_peturb@meta.data %>%
  dplyr::select(true_sgRNA, all_of(score_cols)) %>%
  tidyr::pivot_longer(cols = all_of(score_cols), names_to = "Module", values_to = "Score")

write.csv(scores_long, file.path(OUT_TABS, "module_score_long.csv"), row.names = FALSE)

# Stats vs NTgRNA1
stats <- data.frame(Module = character(), KO = character(), n_KO = integer(), n_NT = integer(),
                    KO_median = numeric(), NT_median = numeric(), p_value = numeric(),
                    stringsAsFactors = FALSE)

for (mod in score_cols) {
  d0 <- PS_peturb@meta.data[, c("true_sgRNA", mod), drop = FALSE]
  colnames(d0) <- c("true_sgRNA", "Score")

  nt <- d0 %>% filter(true_sgRNA == "NTgRNA1") %>% pull(Score)
  nt_med <- median(nt, na.rm = TRUE)

  for (ko in setdiff(unique(d0$true_sgRNA), "NTgRNA1")) {
    kk <- d0 %>% filter(true_sgRNA == ko) %>% pull(Score)
    if (length(nt) < 5 || length(kk) < 5) next

    p <- suppressWarnings(wilcox.test(kk, nt)$p.value)
    stats <- rbind(stats, data.frame(
      Module = mod,
      KO = ko,
      n_KO = length(kk),
      n_NT = length(nt),
      KO_median = median(kk, na.rm = TRUE),
      NT_median = nt_med,
      p_value = p,
      stringsAsFactors = FALSE
    ))
  }
}

stats$p_adj <- p.adjust(stats$p_value, method = "BH")
write.csv(stats, file.path(OUT_TABS, "module_score_stats_vs_NT.csv"), row.names = FALSE)

# -----------------------------
# 9) Plot: module score distributions by KO (violin)
# -----------------------------
pdf(file.path(OUT_PLOTS, "module_scores_violin_by_KO.pdf"), width = 14, height = 6)
for (mod in score_cols) {
  p <- VlnPlot(PS_peturb, features = mod, group.by = "true_sgRNA", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = mod, x = "KO (true_sgRNA)")
  print(p)
}
dev.off()

message("Done. Key outputs written to: ", file.path(PROJECT_ROOT, "results"))
