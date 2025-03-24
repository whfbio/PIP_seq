# ------------------------------------------------------------------------------
# Title: Figure 6 - Summary of Perturbation Results Across Modules
# Author: Haofei Wang
# Date: March 24, 2025
# Description:
# This script visualizes the effect of perturbations using dot plots and module
# summary scores. It uses fold-change and p-values to highlight significantly 
# affected genes and modules.
#
# Dependencies:
# - ggplot2
# - dplyr
# - ggrepel
# - reshape2
#
# Usage:
# Load the preprocessed perturbation result CSV files before execution.
# ------------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(reshape2)

source("./scripts/functions.R")
source("./PS_16/Reference_data_analysis/fig_colors.R")

###Figure6F
generank<-read.csv("./rObjs/gene_rank.csv", header = F)
result.sum<-read.csv("./plots/Perturb_result_summary.csv", row.names = 1)
results.filter <- result.sum %>%
  filter(Measurement %in% c("vCM1","vCM4","GO:0055023","GO:0055017","GO:0061337", "GO:0086042", "GO:1903779", "GO:1990584", "ActivationEarlyheart contraction", "ActivationMidregulation of membrane potential","Channel and handling","Mitochondria","Myofilament"))
results.filter$Gene <- factor(results.filter$Gene, levels = generank$V1)

module_rank<-c("vCM1","vCM4","GO:0055023","GO:0055017","GO:0061337", "GO:0086042", "GO:1903779", "GO:1990584", "ActivationEarlyheart contraction", "ActivationMidregulation of membrane potential","Channel and handling","Mitochondria","Myofilament")
# Convert the 'Measurement' column to a factor with the specified levels
results.filter$Measurement <- factor(results.filter$Measurement, levels = module_rank)

pdf("./plots/PIP_dotplot_filter.pdf", width = 16, height = 6)

ggplot(results.filter, aes(x = Gene, y = Measurement, fill = NormalizedCappedFoldChange, size = abs(NormalizedFoldChange))) +
  geom_point(shape = 21,  # This shape allows for both fill and border color
             aes(color = ifelse(P_Value < 0.05, "red", NA)),  # Conditional border color
             stroke = 1) +  # Border thickness
  scale_fill_gradient2(low = "#568478", high = "#ff8d82") +  # Adjust the color gradient for fill
  scale_color_manual(values = c("black", NA)) +  # Define manual color scale for borders
  scale_size(range = c(1, 6)) +  # Adjust dot sizes
  theme_minimal() +
  labs(fill = "Normalized Fold Change", size = "Normalized Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt x labels for better visibility
        legend.position = "right",
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())

dev.off()


###Figure 6G
###produce perturbed module score
result.sum<-read.csv("./plots/Perturb_result_summary.csv", row.names = 1)
df<-result.sum[result.sum$Measurement%in% c(Myofibril,Electrophysiology,Metabolism),]
df.sig<-result.sum[result.sum$P_Value <=0.05 & result.sum$FoldChange < 0 ,]
result.df<-data.frame(Gene = as.character(), Module=as.numeric())
for (target in unique(df$Gene)){
  df.temp<-df.sig[df.sig$Gene==target,]
  temp.df<-data.frame(Gene = target, Module=nrow(df.temp))
  result.df<-rbind(result.df,temp.df)
}
# Ensure the Module column is treated as a factor
result.df$Module <- as.factor(result.df$Module)
color.schme<-c("0"="#292f56","1"="#00a3a4","2"="#00a3a4","3"="#00bca1","4"="#00bca1","5"="#00d493","6"="#00d493","7"="#69e882","8"="#69e882","9"="#acfa70","10"="#acfa70")
library(ggplot2)
library(ggrepel)
ggplot(result.df, aes(x = Gene, y = Module)) +
  geom_point(color = 'blue') +
  geom_text_repel(aes(label = Gene), 
                  hjust = 1, vjust = 1) +
  scale_color_manual(color.schme) +  # Gradient for non-zero values
  ggtitle('Scatter Plot of Gene Values') +
  xlab('Gene') +
  ylab('Value') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "white"),  # Ensure axis ticks are white
    axis.line = element_line(color = "white"),   # Ensure axis lines are white
    plot.background = element_rect(fill = "black"),  # Set plot background to black
    panel.background = element_rect(fill = "black"), # Set panel background to black
    axis.text = element_text(color = "white"),       # Axis text in white
    axis.title = element_text(color = "white"),      # Axis titles in white
    plot.title = element_text(color = "white")       # Plot title in white
  )


# Create a vector of viridis colors for Module scores from 0 to 10
viridis_colors <- c(
  "0" = "#440154",  # Equivalent of viridis(0)
  "1" = "#482878",
  "2" = "#3E4A89",
  "3" = "#31688E",
  "4" = "#26828E",
  "5" = "#1F9E89",
  "6" = "#35B779",
  "7" = "#6CCE59",
  "8" = "#B4DD2C",
  "9" = "#FDE725",
  "10" = "#FDE725" # Repeat the same for the highest score
)

# Ensure the Module column is treated as a factor
result.df$Module <- as.factor(result.df$Module)

# Plot using the viridis color list for the points
pdf("../plots/Screen_summary_label.pdf", width = 8, height = 5)
ggplot(result.df, aes(x = Gene, y = Module)) +
  geom_point(aes(color = Module), size = 4) +  # Color the points based on Module score
  #geom_text_repel(aes(label = Gene), 
  #                hjust = 1, vjust = 1, 
  #                color = "white", 
  #                family = "Arial",     # Set font to Arial
  #                fontface = "italic") +  # Italicize the labels # Keep the labels in white
  scale_color_manual(values = viridis_colors) +  # Use manual color mapping for Module scores
  ggtitle('Scatter Plot of Gene Values') +
  xlab('Gene') +
  ylab('Value') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "white"),  # Ensure axis ticks are white
    axis.line = element_line(color = "white"),   # Ensure axis lines are white
    plot.background = element_rect(fill = "black"),  # Set plot background to black
    panel.background = element_rect(fill = "black"), # Set panel background to black
    axis.text = element_text(color = "white"),       # Axis text in white
    axis.title = element_text(color = "white"),      # Axis titles in white
    plot.title = element_text(color = "white")       # Plot title in white
  )
dev.off()

###Figure6H
###distance plot for each sgRNA
# install.packages("tidyr")  # if not installed
library(tidyr)
library(dplyr)  # for the pipe (%>%)
result.sum<-read.csv("./plots/Perturb_result_summary.csv", row.names = 1)
result.slim<-result.sum[, c("Gene","Measurement","FoldChange")]
# Convert from long to wide
df_wide <- result.slim %>%
  pivot_wider(
    names_from  = Measurement,  # columns come from go_category
    values_from = FoldChange         # cell values come from value
  )
df_wide<-as.data.frame(df_wide)
rownames(df_wide) <- df_wide$Gene
df_wide<-df_wide[,-1]
df_wide_slim<-df_wide[, module_rank]
df_wide_slim<-df_wide_slim[-10, ]
df_wide_slim<-df_wide_slim[result.df[result.df$color_category == "non-zero", ]$Gene,]
df_wide_slim<-df_wide_slim[-8,]
# View the resulting data frame
print(df_wide)


# Compute the distance matrix using Euclidean distance.
distance_matrix <- dist(df_wide_slim, method = "euclidean")

# Perform hierarchical clustering on the distance matrix.
hc <- hclust(distance_matrix, method = "complete")  # you can choose other methods like "average", "single", etc.

# Plot the dendrogram.
plot(hc, main = "Dendrogram of Genes", xlab = "Genes", sub = "", cex = 0.9)



# Convert the distance object to a matrix for heatmap plotting
d_mat <- as.matrix(distance_matrix)

library(viridis)
library(pheatmap)
pdf("../../../../Manuscript/Figure 6/sgRNA_distance_candidate_viridis.pdf")
pheatmap(d_mat, 
         main = "Gene Distance Heatmap",
         color = rev(viridis(50)))
dev.off()