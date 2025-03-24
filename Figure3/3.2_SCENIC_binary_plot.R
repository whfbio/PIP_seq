# ------------------------------------------------------------------------------
# Title: SCENIC Binary Regulon Activity Visualization
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script loads and processes SCENIC binary regulon activity for Xenium data,
# generates heatmaps of transcription factor (TF) activity, and annotates them
# with cell type and timepoint metadata.
#
# Key Functions:
# - Loads SCENIC binary regulon activity
# - Maps regulon activity to Xenium cell identities
# - Creates heatmaps to visualize regulon activity across cell types and timepoints
#
# Dependencies:
# - Seurat, ggplot2, dplyr, ComplexHeatmap, circlize, reshape2
#
# Output:
# - Regulon activity matrix saved as CSV
# - Binary regulon heatmaps exported as PDFs
#
# Usage:
# - Run this script in an R environment where SCENIC and Xenium objects are available.
# ------------------------------------------------------------------------------

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(readr)

# ------------------------------------------------------------------------------
# Step 1: Load Processed SCENIC and Xenium Data
# ------------------------------------------------------------------------------
ref.scenic <- readRDS("./ref.scenic.xenium.rds")
xenium <- readRDS("./xenium.all.0828.rds")

# Load SCENIC binary regulon activity
scenicOptions <- readRDS("./int/scenicOptions.Rds")
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
binaryRegulonActivity <- binaryRegulonActivity[onlyNonDuplicatedExtended(rownames(binaryRegulonActivity)),]

# Load Xenium metadata mapping
xenium_snRNA_mapping <- read_csv("xenium_snRNA_mapping.csv")

# ------------------------------------------------------------------------------
# Step 2: Map SCENIC Binary Regulon Activity to Xenium Cells
# ------------------------------------------------------------------------------
# Transform binary regulon activity
binaryRegulonActivity_t <- as.data.frame(t(binaryRegulonActivity))
binaryRegulonActivity_t$SCENIC <- rownames(binaryRegulonActivity_t)

# Merge with Xenium cell IDs
mapped_data <- dplyr::right_join(binaryRegulonActivity_t, xenium_snRNA_mapping, by = "SCENIC")
rownames(mapped_data) <- mapped_data$rownames_xenium

# Clean up columns
mapped_data <- mapped_data %>% select(-SCENIC, -`xenium$OriginalCID_snRNA`, -rownames_xenium)

# Convert to matrix
binaryRegulonActivity_x <- as.matrix(t(mapped_data))

# ------------------------------------------------------------------------------
# Step 3: Compute Regulon Activity by Cell Type
# ------------------------------------------------------------------------------
cellInfo <- data.frame(xenium@meta.data)
combinedFactor <- factor(paste(cellInfo$cytoSPACE.final.anno, cellInfo$time_point, sep="_"))

# Compute mean binary regulon activity per cell type
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), combinedFactor), function(cells) {
  if (length(cells) > 1) {
    return(rowMeans(binaryRegulonActivity_x[, cells]))
  } else if (length(cells) == 1) {
    return(binaryRegulonActivity_x[, cells])
  } else {
    return(NA)
  }
})

# Save processed regulon activity
write_csv(as.data.frame(regulonActivity_byCellType), "./regulonActivity_byCellType_xenium.csv")

# ------------------------------------------------------------------------------
# Step 4: Heatmap Visualization
# ------------------------------------------------------------------------------
# Define color scale
my_color_palette <- colorRampPalette(c("white", "#4e8872"))(100)
color_mapping <- colorRamp2(seq(0, 1, length.out = length(my_color_palette)), my_color_palette)

# Convert regulon activity data to ordered format
data_ordered <- regulonActivity_byCellType[, order(colnames(regulonActivity_byCellType))]

# Define heatmap annotations
timepoint_colors <- c("_p0" = "#FDF6B5", "_p7" = "#F9C480", "_p14" = "#F28D6F", "_p21" = "#E24C80")
colnames_split <- strsplit(colnames(data_ordered), "_")
cell_types <- sapply(colnames_split, function(x) x[1])
time_points <- factor(sapply(colnames_split, function(x) paste0("_", x[2])), levels = names(timepoint_colors))

# Heatmap annotations
annotations_combined <- HeatmapAnnotation(
  CellType = cell_types,
  TimePoint = time_points,
  col = list(TimePoint = timepoint_colors),
  annotation_legend_param = list(CellType = list(title = "Cell Type"), TimePoint = list(title = "Time Point"))
)

# Generate heatmap
heatmap_object <- Heatmap(
  data_ordered, 
  name = "Regulon activity", 
  col = colorRamp2(c(min(data_ordered), max(data_ordered)), c("white", "#4699B7")),
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_columns = FALSE, 
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4),
  top_annotation = annotations_combined
)

# Save heatmap as PDF
pdf("./plot/binarized_regulon_activities_Xenium.pdf")
draw(heatmap_object)
dev.off()


