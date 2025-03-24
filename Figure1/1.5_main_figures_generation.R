# ------------------------------------------------------------------------------
# Title: Figure 1 - Cell Count and Cell Composition Visualizations
# Author: Haofei Wang
# Date: March 24, 2025
# Description:
# This script generates plots for Figure 1 of the perturb-seq study, including:
# - Cell count barplots across timepoints (Xenium and snRNA-seq)
# - Doughnut plots of cell composition for major cell types (lv1, Fib, CM, EC)
#
# Dependencies:
# - Seurat
# - ggplot2
# - patchwork
# - pheatmap
# - cowplot
#
# Usage:
# Run in a working directory containing the necessary `ref` and `ref.xenium` objects.
# ------------------------------------------------------------------------------


library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(gridExtra)
library(cowplot)
library(grid)
library(pheatmap)
library(ggrepel)

source("./PS_16/Reference_data_analysis/fig_colors.R")
source("./PS_38/analysis/scripts/functions.R")

ref<-readRDS("./robjs/ref_man_final.rds")
ref.xenium<-readRDS("./ref.xenium.with12cellcharter.rds")

###Figure 1D
###Barplot of Cell Number for Xenium and snRNA-seq cell number
df <- table(ref$time_point) %>% as.data.frame()
df$Var1 <- factor(df$Var1, levels = c("p0", "p7", "p14","p21"))
df$Var2<-c(46747,98188,114524,143709)
# Melt the data to long format
df_melted <- melt(df, id.vars = "Var1", variable.name = "Variable", value.name = "Value")

# Create a new column to control where the pattern should apply (for example only for Var2)
df_melted$Pattern_Flag <- ifelse(df_melted$Variable == "Var2", "pattern", "solid")
pdf("../Figure 1/cell_number_overall.pdf")
# Create the bar plot with patterns and fill color based on Var1
ggplot(df_melted, aes(x = Var1, y = Value, fill = Var1, pattern = Pattern_Flag)) +
  geom_bar_pattern(stat = "identity", position = "dodge", 
                   color = "black", 
                   pattern_density = 0.1, 
                   pattern_spacing = 0.05, 
                   pattern_fill = "white", 
                   pattern_colour = "black") +
  scale_fill_manual(values = time_pal) +  # Color based on Var1
  scale_pattern_manual(values = c("solid" = NA, "pattern" = "circle")) +  # Pattern for Var2
  labs(title = "Barplot with Color Based on Var1 and Pattern for Var2", 
       x = "Category", 
       y = "Values") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

###Figure 1I
##plot for overall cell lv1 percentage
# Create a dataframe with random proportions
data <- data.frame(Time = ref.xenium$time_point, Celltype = ref.xenium$lv1_annotation) %>%
  group_by(Time, Celltype) %>%
  summarise(Count = n()) %>%  # Count the number of occurrences of each cell type
  # Calculate the proportion of each cell type within each time group
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  ungroup()
# Data for the inner pie chart
p0_df <- filter(data, data$Time =="p0")
p7_df <- filter(data, data$Time =="p7")
p14_df <- filter(data, data$Time =="p14")
p21_df <- filter(data, data$Time =="p21")
###use width 0.15 for CM, 0.3 for endo
pdf("../../Figure 1/xenium_cellprop_lv1_doughunt.pdf")
# Create the stacked pie chart
ggplot() +
  # Inner pie chart
  geom_bar(data = p0_df, aes(x = 1/sum(p0_df$Count)*sum(p0_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p7_df, aes(x = 1/sum(p0_df$Count)*sum(p7_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p14_df, aes(x = 1/sum(p0_df$Count)*sum(p14_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p21_df, aes(x = 1/sum(p0_df$Count)*sum(p21_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  # Convert to polar coordinates to make it circular
  coord_polar(theta = "y") +
  scale_fill_manual(values = xenium_lv1_pal) +
  # Set limits for the x-axis to control the spacing between the layers
  xlim(0.5, (1/sum(p0_df$Count)*sum(p21_df$Count)+0.5)) +
  theme_void() +  # Remove background and axes
  theme(legend.position = "right")  # Optional: Position the legend
dev.off()


###Figure 1K
ref.xenium.subset<-subset(ref.xenium, lv1_annotation == "Fib")
# Create a dataframe with random proportions
data <- data.frame(Time = ref.xenium.subset$time_point, Celltype = ref.xenium.subset$cytoSPACE.final.anno) %>%
  group_by(Time, Celltype) %>%
  summarise(Count = n()) %>%  # Count the number of occurrences of each cell type
  # Calculate the proportion of each cell type within each time group
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  ungroup()
# Data for the inner pie chart
p0_df <- filter(data, data$Time =="p0")
p7_df <- filter(data, data$Time =="p7")
p14_df <- filter(data, data$Time =="p14")
p21_df <- filter(data, data$Time =="p21")
pdf("../Figure 1/xenium_cellprop_fib_doughunt.pdf")
# Create the stacked pie chart
ggplot() +
  # Inner pie chart
  geom_bar(data = p0_df, aes(x = 1/sum(p0_df$Count)*sum(p0_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p7_df, aes(x = 1/sum(p0_df$Count)*sum(p7_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p14_df, aes(x = 1/sum(p0_df$Count)*sum(p14_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p21_df, aes(x = 1/sum(p0_df$Count)*sum(p21_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  # Convert to polar coordinates to make it circular
  coord_polar(theta = "y") +
  scale_fill_manual(values = xenium_fib_lv2_pal) +
  # Set limits for the x-axis to control the spacing between the layers
  xlim(0.5, (1/sum(p0_df$Count)*sum(p21_df$Count)+0.5)) +
  theme_void() +  # Remove background and axes
  theme(legend.position = "right")  # Optional: Position the legend
dev.off()


ref.xenium.subset<-subset(ref.xenium, lv1_annotation == "CM")
# Create a dataframe with random proportions
data <- data.frame(Time = ref.xenium.subset$time_point, Celltype = ref.xenium.subset$cytoSPACE.final.anno) %>%
  group_by(Time, Celltype) %>%
  summarise(Count = n()) %>%  # Count the number of occurrences of each cell type
  # Calculate the proportion of each cell type within each time group
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  ungroup()
# Data for the inner pie chart
p0_df <- filter(data, data$Time =="p0")
p7_df <- filter(data, data$Time =="p7")
p14_df <- filter(data, data$Time =="p14")
p21_df <- filter(data, data$Time =="p21")
pdf("../Figure 1/xenium_cellprop_CM_doughunt.pdf")
# Create the stacked pie chart
ggplot() +
  # Inner pie chart
  geom_bar(data = p0_df, aes(x = 1/sum(p0_df$Count)*sum(p0_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.15) +
  geom_bar(data = p7_df, aes(x = 1/sum(p0_df$Count)*sum(p7_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.15) +
  geom_bar(data = p14_df, aes(x = 1/sum(p0_df$Count)*sum(p14_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.15) +
  geom_bar(data = p21_df, aes(x = 1/sum(p0_df$Count)*sum(p21_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.15) +
  # Convert to polar coordinates to make it circular
  coord_polar(theta = "y") +
  scale_fill_manual(values = xenium_cm_lv2_pal) +
  # Set limits for the x-axis to control the spacing between the layers
  xlim(0.5, (1/sum(p0_df$Count)*sum(p21_df$Count)+0.5)) +
  theme_void() +  # Remove background and axes
  theme(legend.position = "right")  # Optional: Position the legend
dev.off()

ref.xenium.subset<-subset(ref.xenium, lv1_annotation == "EC")
# Create a dataframe with random proportions
data <- data.frame(Time = ref.xenium.subset$time_point, Celltype = ref.xenium.subset$cytoSPACE.final.anno) %>%
  group_by(Time, Celltype) %>%
  summarise(Count = n()) %>%  # Count the number of occurrences of each cell type
  # Calculate the proportion of each cell type within each time group
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  ungroup()
# Data for the inner pie chart
p0_df <- filter(data, data$Time =="p0")
p7_df <- filter(data, data$Time =="p7")
p14_df <- filter(data, data$Time =="p14")
p21_df <- filter(data, data$Time =="p21")
pdf("../Figure 1/xenium_cellprop_EC_doughunt.pdf")
# Create the stacked pie chart
ggplot() +
  # Inner pie chart
  geom_bar(data = p0_df, aes(x = 1/sum(p0_df$Count)*sum(p0_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p7_df, aes(x = 1/sum(p0_df$Count)*sum(p7_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p14_df, aes(x = 1/sum(p0_df$Count)*sum(p14_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p21_df, aes(x = 1/sum(p0_df$Count)*sum(p21_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  # Convert to polar coordinates to make it circular
  coord_polar(theta = "y") +
  scale_fill_manual(values = xenium_EC_lv2_pal) +
  # Set limits for the x-axis to control the spacing between the layers
  xlim(0.5, (1/sum(p0_df$Count)*sum(p21_df$Count)+0.5)) +
  theme_void() +  # Remove background and axes
  theme(legend.position = "right")  # Optional: Position the legend
dev.off()



####plot for sub celltype percentage
ref.xenium.subset<-subset(ref.xenium, lv1_annotation %in% c("Blood","SMC","Epic","Neur","Peri"))
# Create a dataframe with random proportions
data <- data.frame(Time = ref.xenium.subset$time_point, Celltype = ref.xenium.subset$cytoSPACE.final.anno) %>%
  group_by(Time, Celltype) %>%
  summarise(Count = n()) %>%  # Count the number of occurrences of each cell type
  # Calculate the proportion of each cell type within each time group
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  ungroup()


##plot for overall cell lv1 percentage
# Create a dataframe with random proportions
data <- data.frame(Time = ref.xenium$time_point, Celltype = ref.xenium$lv1_annotation) %>%
  group_by(Time, Celltype) %>%
  summarise(Count = n()) %>%  # Count the number of occurrences of each cell type
  # Calculate the proportion of each cell type within each time group
  group_by(Time) %>%
  mutate(Proportion = Count / sum(Count) * 100) %>%
  ungroup()
# Data for the inner pie chart
p0_df <- filter(data, data$Time =="p0")
p7_df <- filter(data, data$Time =="p7")
p14_df <- filter(data, data$Time =="p14")
p21_df <- filter(data, data$Time =="p21")
###use width 0.15 for CM, 0.3 for endo
pdf("../../Figure 1/xenium_cellprop_lv1_doughunt.pdf")
# Create the stacked pie chart
ggplot() +
  # Inner pie chart
  geom_bar(data = p0_df, aes(x = 1/sum(p0_df$Count)*sum(p0_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p7_df, aes(x = 1/sum(p0_df$Count)*sum(p7_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p14_df, aes(x = 1/sum(p0_df$Count)*sum(p14_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  geom_bar(data = p21_df, aes(x = 1/sum(p0_df$Count)*sum(p21_df$Count), y = Proportion, fill = Celltype), 
           stat = "identity", width = 0.3) +
  # Convert to polar coordinates to make it circular
  coord_polar(theta = "y") +
  scale_fill_manual(values = xenium_lv1_pal) +
  # Set limits for the x-axis to control the spacing between the layers
  xlim(0.5, (1/sum(p0_df$Count)*sum(p21_df$Count)+0.5)) +
  theme_void() +  # Remove background and axes
  theme(legend.position = "right")  # Optional: Position the legend
dev.off()

###Figure 1L
###scatter plot showing the trend of proliferative cells in the data
meta<-ref.xenium@meta.data
result <- data.frame(
  "celltype" = character(),
  "p0" = numeric(),
  "p7" = numeric(),
  "p14" = numeric(),
  "p21" = numeric()
)
temp <- data.frame(
  "celltype" = "p.EC",
  "p0" = nrow(meta[meta$cytoSPACE.final.anno == "p.EC" & meta$time_point == "p0",])/nrow(meta[meta$time_point == "p0",])*100,
  "p7" = nrow(meta[meta$cytoSPACE.final.anno == "p.EC" & meta$time_point == "p7",])/nrow(meta[meta$time_point == "p7",])*100,
  "p14" = nrow(meta[meta$cytoSPACE.final.anno == "p.EC" & meta$time_point == "p14",])/nrow(meta[meta$time_point == "p14",])*100,
  "p21" = nrow(meta[meta$cytoSPACE.final.anno == "p.EC" & meta$time_point == "p21",])/nrow(meta[meta$time_point == "p21",])*100
)
result<-rbind(result, temp)
temp <- data.frame(
  "celltype" = "p.Fib",
  "p0" = nrow(meta[meta$cytoSPACE.final.anno == "p.Fib" & meta$time_point == "p0",])/nrow(meta[meta$time_point == "p0",])*100,
  "p7" = nrow(meta[meta$cytoSPACE.final.anno == "p.Fib" & meta$time_point == "p7",])/nrow(meta[meta$time_point == "p7",])*100,
  "p14" = nrow(meta[meta$cytoSPACE.final.anno == "p.Fib" & meta$time_point == "p14",])/nrow(meta[meta$time_point == "p14",])*100,
  "p21" = nrow(meta[meta$cytoSPACE.final.anno == "p.Fib" & meta$time_point == "p21",])/nrow(meta[meta$time_point == "p21",])*100
)
result<-rbind(result, temp)
temp <- data.frame(
  "celltype" = "p.CM",
  "p0" = nrow(meta[meta$cytoSPACE.final.anno == "p.CM" & meta$time_point == "p0",])/nrow(meta[meta$time_point == "p0",])*100,
  "p7" = nrow(meta[meta$cytoSPACE.final.anno == "p.CM" & meta$time_point == "p7",])/nrow(meta[meta$time_point == "p7",])*100,
  "p14" = nrow(meta[meta$cytoSPACE.final.anno == "p.CM" & meta$time_point == "p14",])/nrow(meta[meta$time_point == "p14",])*100,
  "p21" = nrow(meta[meta$cytoSPACE.final.anno == "p.CM" & meta$time_point == "p21",])/nrow(meta[meta$time_point == "p21",])*100
)
result<-rbind(result, temp)
df_long <- pivot_longer(result, 
                        cols = -celltype, 
                        names_to = "timepoint", 
                        values_to = "percentage")

# Optionally, if the timepoint order is not correct, set the factor levels explicitly
df_long$timepoint <- factor(df_long$timepoint, levels = c("p0", "p7", "p14", "p21"))
df_long$time_numeric <- as.numeric(factor(df_long$timepoint, levels = c("p0", "p7", "p14", "p21")))


# Create a line plot with one line per cell type
# Now plot using ggplot2:
celltype_pal<-c("p.CM" = "#CFC8DA","p.EC"="#CEE9F4","p.Fib"="#ffcccb")
pdf("../../../Manuscript/Figure 1/proliferative_cell_num.pdf")
ggplot(df_long, aes(x = time_numeric, y = percentage, group = celltype, fill = celltype)) +
  geom_area(alpha = 0.3, position = "identity") +  # fills the area under each line
  geom_line(aes(color = celltype), size = 1) +       # adds the lines
  geom_point(aes(color = celltype), size = 3) +      # adds the points
  scale_x_continuous(breaks = unique(df_long$time_numeric),
                     labels = unique(df_long$timepoint)) +  # set x-axis labels as T1, T2, etc.
  scale_color_manual(values = celltype_pal) +       # set custom line/point colors
  scale_fill_manual(values = celltype_pal) +   
  labs(title = "Percentage of Cells Across Timepoints by Cell Type",
       x = "Timepoint",
       y = "Percentage of Cells") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # remove major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    axis.line = element_line(color = "black")  # add axis lines
  )
dev.off()

#Figure 1H and 1J
###plot xenium plot for cell types
###plot cell annotation plot
#lv1 for all celltype in heart
ref.xenium<-readRDS("./robjs/ref.xenium.withcellcharter.rds")
timepoint<-"p0"
xenium.obj<-subset(ref.xenium,time_point==timepoint)
fov.list<-list("p0"="fov","p7"="fov.2","p14"="fov.3","p21"="fov.4")
cluster_pal<-cell_lv1.xenium.pal
cluster_pal[]<-"black"
cluster_pal<-c(cluster_pal,"others"="#eeeeee")
p<- list()
# Loop through each time point and create plots
for (celltype in unique(xenium.obj$lv1_annotation)) {
  xenium.temp <- create_custom_cluster_metadata(xenium.obj, "lv1_annotation", celltype)
  
  # Generate the UMAP plot for the current cluster
  p[["all"]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "lv1_annotation", cols = cell_lv1.xenium.pal,
                            boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
    ggtitle(paste0(timepoint,"_all"))
  p[[celltype]] <- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "new_clusters", 
                                boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F, cols = cluster_pal) +
    ggtitle(paste0(timepoint,"_",celltype))
  print(paste0("finish plotting source ", celltype))
  print(paste0("progress ", round((length(p) / length(unique(xenium.temp$lv1_annotation)) * 100)), "%"))
 
}
p_rotate <- lapply(p, rotate_image, rot_angle=90)
pdf(paste0("./plots/",timepoint,"_xenium_celltype_lv1.pdf"),onefile = T, width = 9, height = ceiling(length(unique(xenium.obj$lv1_annotation))/3)*3)
grid.arrange(grobs = p_rotate, ncol = 3)  # Arrange in a 2-column grid
dev.off()


###plot subcluster spatial location
xenium_cm_lv2_pal<-c("Ddc.CM"="#651251", "v.CM"= '#1E0012', "Ankrd1.CM"='#D29ACA', "p.CM" = "#CFC8DA","Slit2.CM"="#B9A2BE" ,"a.CM" ="#9F7498","av.CM"="#72476D")
xenium_ec_lv2_pal<-c("cap.EC" ="#183642","arteriole.EC"="#6DB2BF","endocardial.EC"="#4699B7", "p.EC"="#CEE9F4","lym.EC"="#94CFC9","cap.EC.High.Myo."="#20566E")
xenium_fib_lv2_pal<-c("Gfpt2.Fib"="#cb4154","Col8a1.Fib"="#d55d6c","Lsamp.Fib"="#e07984","Postn.Fib"="#ea949b", "valve.Fib"="#f5b0b3", "p.Fib"="#ffcccb")
xenium_other_lv2_pal <- c("Blood"="#f4895f", "Epic"="#95cf92", "SMC"="#ffd92f", "Neur"="#AAF6D2", "p.Peri"="#b58188","Peri"="#6f1926")
xenium_lv2_pal<-c(xenium_cm_lv2_pal,xenium_ec_lv2_pal,xenium_fib_lv2_pal,xenium_other_lv2_pal)

subcluster<-"Fib"
timepoint <-"p7"
cluster_pal<-c(xenium_fib_lv2_pal,"others"="#eeeeee")
xenium.obj<-subset(ref.xenium,time_point==timepoint&lv1_annotation==subcluster)
p<- list()
# Loop through each time point and create plots
for (celltype in unique(xenium.obj$cytoSPACE.final.anno)) {
  xenium.temp <- create_custom_cluster_metadata(xenium.obj, "cytoSPACE.final.anno", celltype)
  
  # Generate the UMAP plot for the current cluster
  p[["all"]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "cytoSPACE.final.anno", cols = cluster_pal,
                            boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
    ggtitle(paste0(timepoint,"_all"))
  p[[celltype]] <- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "new_clusters", 
                                boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F, cols = cluster_pal) +
    ggtitle(paste0(timepoint,"_",celltype))
  print(paste0("finish plotting source ", celltype))
  print(paste0("progress ", round((length(p) / length(unique(xenium.temp$cytoSPACE.final.anno)) * 100)), "%"))
  
}
p_rotate <- lapply(p, rotate_image, rot_angle=180)
pdf(paste0("./plots/",timepoint,"_xenium_",subcluster,".pdf"),onefile = T, width = length(unique(xenium.obj$cytoSPACE.final.anno)) *4, height = 5)
grid.arrange(grobs = p_rotate, nrow = 1)  # Arrange in a 2-column grid
dev.off()

subcluster<-"CM"
timepoint <-"p7"
cluster_pal<-c(xenium_cm_lv2_pal,"others"="#eeeeee")
xenium.obj<-subset(ref.xenium,time_point==timepoint&lv1_annotation==subcluster)
p<- list()
# Loop through each time point and create plots
for (celltype in unique(xenium.obj$cytoSPACE.final.anno)) {
  xenium.temp <- create_custom_cluster_metadata(xenium.obj, "cytoSPACE.final.anno", celltype)
  
  # Generate the UMAP plot for the current cluster
  p[["all"]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "cytoSPACE.final.anno", cols = cluster_pal,
                            boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
    ggtitle(paste0(timepoint,"_all"))
  p[[celltype]] <- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "new_clusters", 
                                boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F, cols = cluster_pal) +
    ggtitle(paste0(timepoint,"_",celltype))
  print(paste0("finish plotting source ", celltype))
  print(paste0("progress ", round((length(p) / length(unique(xenium.temp$cytoSPACE.final.anno)) * 100)), "%"))
  
}
p_rotate <- lapply(p, rotate_image, rot_angle=180)
pdf(paste0("./plots/",timepoint,"_xenium_",subcluster,".pdf"),onefile = T, width = length(unique(xenium.obj$cytoSPACE.final.anno)) *4, height = 5)
grid.arrange(grobs = p_rotate, nrow = 1)  # Arrange in a 2-column grid
dev.off()

subcluster<-"EC"
timepoint <-"p7"
cluster_pal<-c(xenium_ec_lv2_pal,"others"="#eeeeee")
xenium.obj<-subset(ref.xenium,time_point==timepoint&lv1_annotation==subcluster)
p<- list()
# Loop through each time point and create plots
for (celltype in unique(xenium.obj$cytoSPACE.final.anno)) {
  xenium.temp <- create_custom_cluster_metadata(xenium.obj, "cytoSPACE.final.anno", celltype)
  
  # Generate the UMAP plot for the current cluster
  p[["all"]]<- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "cytoSPACE.final.anno", cols = cluster_pal,
                            boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F) +
    ggtitle(paste0(timepoint,"_all"))
  p[[celltype]] <- ImageDimPlot(xenium.temp, fov = fov.list[[timepoint]], nmols = 20000, group.by = "new_clusters", 
                                boundaries = "segmentation", border.color = NA, dark.background = F, crop = T, axes = F, cols = cluster_pal) +
    ggtitle(paste0(timepoint,"_",celltype))
  print(paste0("finish plotting source ", celltype))
  print(paste0("progress ", round((length(p) / length(unique(xenium.temp$cytoSPACE.final.anno)) * 100)), "%"))
  
}
p_rotate <- lapply(p, rotate_image, rot_angle=180)
pdf(paste0("./plots/",timepoint,"_xenium_",subcluster,".pdf"),onefile = T, width = length(unique(xenium.obj$cytoSPACE.final.anno)) *4, height = 5)
grid.arrange(grobs = p_rotate, nrow = 1)  # Arrange in a 2-column grid
dev.off()


