# ------------------------------------------------------------------------------
# Title: SCENIC Integration with Xenium Data
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script integrates SCENIC results with Xenium spatial transcriptomics data,
# processes transcription factor activity, and performs hierarchical clustering.
#
# Key Functions:
# - Maps single-cell RNA-seq (scRNA-seq) data to spatial transcriptomics (ST)
# - Computes SCENIC regulon activity
# - Generates heatmaps and hierarchical clustering analysis for TF activity
#
# Dependencies:
# - Seurat, ggplot2, dplyr, ComplexHeatmap, AUCell, circlize, RColorBrewer
#
# Output:
# - SCENIC-mapped Xenium Seurat object
# - Regulon activity clustering
# - Hierarchical heatmap of regulon activity
#
# Usage:
# - Run this script in an R environment where the required Seurat and SCENIC
#   objects are accessible.
# ------------------------------------------------------------------------------

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(AUCell)
library(circlize)
library(RColorBrewer)
library(reshape2)
library(Scillus)
library(SCENIC)
library(arrow)

# ------------------------------------------------------------------------------
# Load and Process SCENIC and Xenium Data
# ------------------------------------------------------------------------------

# Load reference SCENIC object and Xenium data
ref.scenic <- readRDS("./ref.scenic.rds")
xenium <- readRDS("./xenium.all.0828.rds")

# Identify common cells between SCENIC and Xenium
common_cells <- intersect(colnames(ref.scenic), xenium$OriginalCID_snRNA)
ref.scenic <- subset(ref.scenic, cells = common_cells)

# Map SCENIC metadata to Xenium
xenium_snRNA_mapping <- data.frame(rownames_xenium = colnames(xenium), 
                                   OriginalCID_snRNA = xenium$OriginalCID_snRNA,
                                   cytoSPACE.final.anno = xenium$cytoSPACE.final.anno)
xenium_snRNA_mapping <- dplyr::left_join(xenium_snRNA_mapping, 
                                         data.frame(rownames_xenium = rownames(xenium@meta.data), 
                                                    cytoSPACE.final.anno = xenium$cytoSPACE.final.anno), 
                                         by = "rownames_xenium")

# Add mapped metadata to SCENIC reference object
ref.scenic <- AddMetaData(ref.scenic, metadata = xenium_snRNA_mapping$rownames_xenium, col.name = "xenium_barcode")
ref.scenic <- AddMetaData(ref.scenic, metadata = xenium_snRNA_mapping$cytoSPACE.final.anno, col.name = "cytoSPACE.final.anno")

# Save updated SCENIC reference object
write_csv(xenium_snRNA_mapping, "xenium_snRNA_mapping.csv")
saveRDS(ref.scenic, "./ref.scenic.xenium.rds")

# Create scenic_mat_mapped_data.feather for h5ad Xenium File integration.
scenic_mat <- ref.scenic@assays$SCENIC@data
xenium_snRNA_mapping <- read_csv("xenium_snRNA_mapping.csv")
dim(scenic_mat)
scenic_mat_T <- as.data.frame(t(scenic_mat))
scenic_mat_T$SCENIC <- rownames(scenic_mat_T)
scenic_mat_mapped_data <- dplyr::right_join(scenic_mat_T, xenium_snRNA_mapping, by = "SCENIC")
rownames(scenic_mat_mapped_data) <- scenic_mat_mapped_data$rownames_xenium
scenic_mat_mapped_data$SCENIC<-NULL
scenic_mat_mapped_data$`xenium$OriginalCID_snRNA`<-NULL
scenic_mat_mapped_data$rownames_xenium<-NULL
scenic_mat_mapped_data <- as.matrix(scenic_mat_mapped_data)
scenic_mat_mapped_data <- t(scenic_mat_mapped_data)
write_feather(as.data.frame(scenic_mat_mapped_data), "scenic_mat_mapped_data.feather")

# ------------------------------------------------------------------------------
# Create Heatmap of SCENIC Regulon Activity
# ------------------------------------------------------------------------------

# Extract SCENIC expression data
scenic_data <- ref.scenic@assays$SCENIC@data

# Reorder column names to match Xenium metadata
ref_to_xenium_barcode <- setNames(xenium_snRNA_mapping$rownames_xenium, xenium_snRNA_mapping$OriginalCID_snRNA)
colnames(scenic_data) <- make.unique(ref_to_xenium_barcode[colnames(scenic_data)])

# Add SCENIC assay to Xenium object
xenium[["SCENIC"]] <- CreateAssay5Object(data = as(scenic_data, "dgCMatrix"))

# Generate hierarchical clustering
hclust_matrix <- scale(t(scenic_data)) %>% t()
gene_hclust <- hclust(dist(hclust_matrix), method = "ward.D2")

# Define heatmap colors
moma_pal <- c("#e35496", "#FFFFFF", "#a3c452")

# Generate heatmap
ht.moma <- plot_heatmap(dataset = xenium, 
                        markers = rownames(scenic_data)[gene_hclust$order],
                        sort_var = c("cytoSPACE.final.anno", "time_point"),
                        anno_var = c("cytoSPACE.final.anno", "time_point"),
                        hm_colors = moma_pal,
                        row_font_size = 4)

# Save heatmap as PDF
pdf("./plot/scenic_heatmap.pdf", width = 20, height = 15)
print(ht.moma)
dev.off()

