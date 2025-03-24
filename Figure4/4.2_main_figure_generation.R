# ------------------------------------------------------------------------------
# Title: Figure 4 - Maturation Feature Mapping and GO Analysis
# Author: Haofei Wang
# Date: March 24, 2025
# Description:
# This script generates plots related to maturation features across Xenium timepoints.
# Includes:
# - ImageFeaturePlot for maturation score
# - Zoomed-in views for selected regions
# - Violin plots for maturation across timepoints
# - Spatial neighbor composition analysis
#
# Dependencies:
# - Seurat
# - ggplot2
# - viridis
# - ComplexHeatmap
# - patchwork
#
# Usage:
# Requires Xenium objects with Maturation_Features1 score computed.
# ------------------------------------------------------------------------------


library(Seurat)
library(clusterProfiler)
library(dplyr)
library(org.Mm.eg.db)
library(ReactomePA)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(patchwork)
library(gridExtra)
library(grid)
library(ComplexHeatmap)
library(tidyverse)
###GO plot for WGCNA analysis
setwd("./PS_16/Reference_data_analysis/WGCNA")
source("./scripts/functions.R")
source("./PS_16/Reference_data_analysis/fig_colors.R")

###Figure 4C
###draw the maturation score on heart
# Generate a color palette
library(viridis)
color_palette <- inferno(20)

# Load Xenium datasets
timepoints <- c("p0", "p7", "p14", "p21")
xenium.obj <- lapply(timepoints, function(tp) readRDS(paste0("./robjs/ref.xenium.", tp, ".1011.rds")))
names(xenium.obj) <- timepoints
cm.type<-list("p0"=c("v.CM","Slit2.CM","p.CM","Ankrd1.CM","Ddc.CM","a.CM"),"p7"=c("v.CM","Slit2.CM","p.CM","Ankrd1.CM","Ddc.CM","a.CM"),"p14"=c("v.CM","Slit2.CM","p.CM","Ankrd1.CM","a.CM"),"p21"=c("v.CM","Slit2.CM","Ankrd1.CM","a.CM"))
heatmap.list<-list()
# Loop over each timepoint
for (timepoint in names(xenium.obj)) {

  xenium.temp<-subset(xenium.obj[[timepoint]], cytoSPACE.final.anno %in% cm.type[[timepoint]])

  # Generate heatmaps only if pathways are found
  heatmap.list[[timepoint]]<-ImageFeaturePlot(xenium.temp, features = "Maturation_Features1", max.cutoff = 10, size = 1, cols = color_palette, border.size = NA)
  
}

# Save all plots to a PDF
pdf(file.path("../../../Manuscript/Figure 4/maturation_temporal.pdf"), onefile = TRUE)
invisible(lapply(heatmap.list, print))  # Use invisible to suppress unnecessary output
dev.off()
##zoom in view
timepoint<-"p7"
xenium.temp<-subset(xenium.obj[[timepoint]], cytoSPACE.final.anno %in% cm.type[[timepoint]])
options(future.globals.maxSize = 1000 * 1024^2) 
cropped.coords <- Crop(xenium.temp[["fov.2"]], y = c(3000, 4000), x = c(500, 2000), coords = "plot")
xenium.temp[["zoom"]] <- cropped.coords
DefaultBoundary(xenium.temp[["zoom"]]) <- "segmentation"
pdf("../../../Manuscript/Figure 4/p7_maturation_zoom.pdf")
ImageFeaturePlot(xenium.temp, fov= "zoom", features = "Maturation_Features1", max.cutoff = 10, size = 1, cols = color_palette, border.size = NA)
dev.off()

###Figure 4D
###violin plot for maturation score
###maturation score in each cm niches in each timepoint

wd<-"./PS_38/analysis/"
setwd(wd)
library(ggpubr)
timepoint<-"p7"
xenium.temp <- readRDS("./robjs/ref.xenium.p0.1011.rds")
Idents(xenium.temp)<-"cytoSPACE.final.anno"
cm_type<-c("v.CM")
xenium.p0.cm <- subset(xenium.temp, idents = cm_type)
xenium.temp <- readRDS("./robjs/ref.xenium.p7.1011.rds")
Idents(xenium.temp)<-"cytoSPACE.final.anno"
cm_type<-c("v.CM")
xenium.p7.cm <- subset(xenium.temp, idents = cm_type)
xenium.temp <- readRDS("./robjs/ref.xenium.p14.1011.rds")
Idents(xenium.temp)<-"cytoSPACE.final.anno"
cm_type<-c("v.CM")
xenium.p14.cm <- subset(xenium.temp, idents = cm_type)
xenium.temp <- readRDS("./robjs/ref.xenium.p21.1011.rds")
Idents(xenium.temp)<-"cytoSPACE.final.anno"
cm_type<-c("v.CM")
xenium.p21.cm <- subset(xenium.temp, idents = cm_type)
# Find the maximum length of the vectors
df1<-data.frame("maturation_score"=xenium.p0.cm@meta.data$Maturation_Features1, "timepoint"="p0")
df2<-data.frame("maturation_score"=xenium.p7.cm@meta.data$Maturation_Features1, "timepoint"="p7")
df3<-data.frame("maturation_score"=xenium.p14.cm@meta.data$Maturation_Features1, "timepoint"="p14")
df4<-data.frame("maturation_score"=xenium.p21.cm@meta.data$Maturation_Features1, "timepoint"="p21")
df<-rbind(df1,df2,df3,df4)
df$timepoint <- factor(df$timepoint, levels = c("p0", "p7", "p14", "p21"))

# Create the violin plot
pdf("../../../Manuscript/Figure 4/maturation_timepoint.pdf", height = 10, width = 6)
ggplot(df, aes(x = timepoint, y = maturation_score, fill = timepoint)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  scale_fill_manual(values = time_pal) +
  stat_compare_means(
    comparisons = list(c("p0", "p7"), c("p7", "p14"), c("p14", "p21")), 
    method = "wilcox.test", 
    label = "p.signif"
  ) +  # Add pairwise comparisons with significance labels
  labs(
    title = "Violin Plot of Maturation Scores with Statistical Tests",
    x = "Timepoint",
    y = "Maturation Score"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
dev.off()

###Figure 4E
neighbour.id <- c("0","1","2","3","4","5","6","7","8","9","10","11")
color_palette <- viridis(length(neighbour.id))
# Create a named vector of colors
element_colors <- setNames(rev(color_palette), c("0","1","2","3","4","5","6","7","8","9","10","11"))
spatial_name<-c("Ddc.CM_spatial_neighbour","cap.EC_spatial_neighbour","v.CM_spatial_neighbour","Peri_spatial_neighbour","cap.EC.High.Myo._spatial_neighbour",
                "Ankrd1.CM_spatial_neighbour","p.CM_spatial_neighbour","Postn.Fib_spatial_neighbour","arteriole.EC_spatial_neighbour",    
                "Slit2.CM_spatial_neighbour","endocardial.EC_spatial_neighbour","p.EC_spatial_neighbour","Blood_spatial_neighbour",
                "Epic_spatial_neighbour","p.Fib_spatial_neighbour","a.CM_spatial_neighbour","valve.Fib_spatial_neighbour","av.CM_spatial_neighbour",           
                "SMC_spatial_neighbour","Neur_spatial_neighbour","Gfpt2.Fib_spatial_neighbour")  
cell.neighbour.p7<-read.csv("../../Experiment/PS_38/analysis/cell distance analysis/xenium_p7_metadata_cellneighbour_v2.csv")
# 0 is the cell type, 1 is the level 1 neighbour, 2 is the level 2 neighbour, 3 is the level 3 neighbour, the rest will be level 4
ref.xenium.p7<-readRDS("../../Experiment/PS_38/analysis/robjs/ref.xenium.p7.1011.rds")
meta.p7<-ref.xenium.p7@meta.data

unique_to_seurat <- setdiff(unique(meta.p7$orig.id), unique(cell.neighbour.p7$orig.id))
if (length(unique_to_seurat)==0){
  scanpy.meta_sorted <- cell.neighbour.p7[match(meta.p7$orig.id, cell.neighbour.p7$orig.id), ]
  for (name in spatial_name){
    ref.xenium.p7<-AddMetaData(ref.xenium.p7, metadata = scanpy.meta_sorted, col.name = name)
  }
} else {
  print("Different index in scanpy metadata and seurat obj!")
}
cropped.coords <- Crop(ref.xenium.p7[["fov.2"]], y = c(3000, 4000), x = c(500, 2000), coords = "plot")
ref.xenium.p7[["zoom"]] <- cropped.coords
DefaultBoundary(ref.xenium.p7[["zoom"]]) <- "segmentation"


p<-list()
for (name in spatial_name){
  p[[name]]<-ImageDimPlot(ref.xenium.p7, fov = "zoom", nmols = 20000, group.by = name, dark.background = F, crop = T, cols = element_colors, boundaries = "segmentation", border.color = NA)
}

p_rotate <- lapply(p, rotate_image, rot_angle=180)
pdf(paste0("../../../Manuscript/Figure 4/p7_zoom_cellneighbour_densityplot_v2.pdf"),onefile = T, width = 7, height = 5)
for (i in 1:length(p_rotate)){
  print(p_rotate[i])
}
dev.off()

# Subset the Seurat object
Idents(ref.xenium.p7) <- "source"
source.type <- c("ventricle_left_compact", "trabecular_right", "ventricle_right_compact", "trabecular_left")
xenium.temp <- subset(ref.xenium.p7, source %in% source.type)
Idents(xenium.temp) <- "cytoSPACE.final.anno"
xenium.temp<- subset(xenium.temp, cytoSPACE.final.anno %in% cm.type[["p7"]])
metadata<-xenium.temp@meta.data

# Initialize the plot list
p <- list()

for (spatial_n in spatial_name) {
  if (!spatial_n %in% colnames(metadata)) {
    warning(paste("Column", spatial_n, "not found in metadata for cell type:", spatial_n))
    next  # Skip to the next iteration if column doesn't exist
  }
  
  # Ensure the grouping column is numeric
  metadata[[spatial_n]] <- as.numeric(as.character(metadata[[spatial_n]]))
  
  # Calculate mean and standard error
  summary_data <- metadata %>%
    group_by(.data[[spatial_n]]) %>%
    summarise(
      mean_expression = mean(Maturation_Features1, na.rm = TRUE),
      sd_expression = sd(Maturation_Features1, na.rm = TRUE),
      n = n()
    ) %>%
    mutate(
      se_expression = sd_expression / sqrt(n),
      spatial_numeric = .data[[spatial_n]]
    )
  
  # Fit a loess model to get smooth curve data
  loess_fit <- loess(mean_expression ~ spatial_numeric, data = summary_data)
  
  # Generate data for the smooth line
  x_vals <- seq(min(summary_data$spatial_numeric), max(summary_data$spatial_numeric), length.out = 100)
  y_vals <- predict(loess_fit, newdata = data.frame(spatial_numeric = x_vals))
  
  # Create a data frame with the smooth line data
  smooth_data <- data.frame(
    x = x_vals,
    y = y_vals
  )
  
  name <- paste0(spatial_n, "_Maturation_score")
  
  # Plot
  p[[name]] <- ggplot(summary_data, aes(x = spatial_numeric, y = mean_expression)) +
    geom_point(color = "blue", alpha = 0.5) +
    geom_line(data = smooth_data, aes(x = x, y = y, color = x), size = 1) +  # Map color to x for gradient
    scale_color_viridis_c(direction = -1) +  # Use viridis color scale
    ylim(1, 3) +  # Set y-axis limits from 1 to 3
    labs(
      title = paste("Maturation Features Across", spatial_n, "Neighborhood Layers"),
      x = "Neighborhood Layer (Distance from Cell Type)",
      y = "Mean Maturation Feature",
      color = "Neighborhood Layer"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),       # Remove grid lines
      axis.line = element_line(color = "black"),  # Add x and y axis lines
      axis.text = element_text(size = 12),        # Adjust axis text size
      axis.title = element_text(size = 14),       # Adjust axis title size
      plot.title = element_text(hjust = 0.5, size = 16)  # Center and size the title
    )
}

pdf(paste0("../Figure 4/p7_cellneighbour_maturation_lineplot_colored_scaled.pdf"),onefile = T, width = 8, height = 5)
for (i in 1:length(p)){
  print(p[[i]])
}
dev.off()

# Initialize the plot list
p <- list()

for (spatial_n in spatial_name) {
  if (!spatial_n %in% colnames(metadata)) {
    warning(paste("Column", spatial_n, "not found in metadata for cell type:", spatial_n))
    next  # Skip to the next iteration if column doesn't exist
  }
  
  # Ensure the grouping column is numeric
  metadata[[spatial_n]] <- as.numeric(as.character(metadata[[spatial_n]]))
  
  # Calculate mean and standard error
  summary_data <- metadata %>%
    group_by(.data[[spatial_n]]) %>%
    summarise(
      mean_expression = mean(Maturation_Features1, na.rm = TRUE),
      sd_expression = sd(Maturation_Features1, na.rm = TRUE),
      n = n()
    ) %>%
    mutate(
      se_expression = sd_expression / sqrt(n),
      spatial_numeric = .data[[spatial_n]]
    )
  
  # Fit a loess model to get smooth curve data
  loess_fit <- loess(mean_expression ~ spatial_numeric, data = summary_data)
  
  # Generate data for the smooth line
  x_vals <- seq(min(summary_data$spatial_numeric), max(summary_data$spatial_numeric), length.out = 100)
  y_vals <- predict(loess_fit, newdata = data.frame(spatial_numeric = x_vals))
  
  # Create a data frame with the smooth line data
  smooth_data <- data.frame(
    x = x_vals,
    y = y_vals
  )
  
  name <- paste0(spatial_n, "_Maturation_score")
  
  # Plot
  p[[name]] <- ggplot(summary_data, aes(x = spatial_numeric, y = mean_expression)) +
    geom_point(color = "blue", alpha = 0.5) +
    geom_line(data = smooth_data, aes(x = x, y = y), color="red",size = 1) +  # Map color to x for gradient
    #scale_color_viridis_c(direction = -1) +  # Use viridis color scale
    ylim(1, 3) +  # Set y-axis limits from 1 to 3
    labs(
      title = paste("Maturation Features Across", spatial_n, "Neighborhood Layers"),
      x = "Neighborhood Layer (Distance from Cell Type)",
      y = "Mean Maturation Feature",
      color = "Neighborhood Layer"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),       # Remove grid lines
      axis.line = element_line(color = "black"),  # Add x and y axis lines
      axis.text = element_text(size = 12),        # Adjust axis text size
      axis.title = element_text(size = 14),       # Adjust axis title size
      plot.title = element_text(hjust = 0.5, size = 16)  # Center and size the title
    )
}

pdf(paste0("../Figure 4/p7_cellneighbour_maturation_lineplot_scaled.pdf"),onefile = T, width = 8, height = 5)
for (i in 1:length(p)){
  print(p[[i]])
}
dev.off()


###Figure 4H
####Ptprm signal visulization
library(viridis)
p7.xenium<-readRDS("./robjs/ref.xenium.p7.direct.1011.rds")
color_palette <- viridis(20)
cropped.coords <- Crop(p7.xenium[["fov.2"]], y = c(3000, 4000), x = c(500, 2000), coords = "plot")
p7.xenium[["zoom"]] <- cropped.coords
DefaultBoundary(p7.xenium[["zoom"]]) <- "segmentation"

ImageFeaturePlot(p7.xenium, fov= "zoom", features = c("commot.user_database.Ptprm.Ptprm", "commot.user_database.PTPRM"), max.cutoff = c(90,90),size = 0.75, cols = color_palette, border.size = NA)
pdf("../../../Manuscript/Figure 4/test.pdf")
ImageFeaturePlot(p7.xenium, fov= "fov.2", features = c("commot.user_database.Ptprm.Ptprm"), max.cutoff = 1,size = 0.75, cols = color_palette, border.size = NA)

ImageDimPlot(p7.xenium, fov= "zoom", group.by = "cytoSPACE.final.anno", cols = cell_color, size = 0.75,border.size = NA)
dev.off()
cell_color<-c("Blood"="#440154","a.CM"="#440154","Postn.Fib"="#440154","p.CM"="#440154","cap.EC"="#f68115","endocardial.EC"="#440154","Gfpt2.Fib"="#440154","Slit2.CM"="#440154",        
              "SMC"="#440154","arteriole.EC"="#440154","Epic"="#440154","av.CM"="#440154","Ankrd1.CM"="#440154","p.Fib"="#440154","cap.EC.High.Myo."="#440154","p.EC"="#440154",
               "valve.Fib"="#440154","Neur"="#440154","Peri"="#440154", "v.CM" ="#51c56a", "Ddc.CM"="#440154", "Col8a1.Fib"="#440154")

