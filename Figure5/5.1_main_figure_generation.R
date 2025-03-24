# ------------------------------------------------------------------------------
# Title: Figure 5 - hdWGCNA Co-expression and Module Hub Analysis
# Author: Haofei Wang
# Date: March 24, 2025
# Description:
# This script performs weighted gene co-expression network analysis (hdWGCNA) on 
# single-cell data, identifies hub genes, calculates eigengenes, and visualizes:
# - Module dendrograms
# - Module connectivity and hub gene plots
# - Module UMAPs and networks
#
# Dependencies:
# - Seurat
# - hdWGCNA
# - WGCNA
# - UCell
# - ggplot2
#
# Usage:
# Load the reference vCM object before running.
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
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
setwd("./PS_16/Reference_data_analysis/WGCNA")
source("./scripts/functions.R")
source("./PS_16/Reference_data_analysis/fig_colors.R")

###Figure 5A-B
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load snRNA-seq dataset
wd<-c("./PS_16/Reference_data_analysis/WGCNA")
seurat_obj <- readRDS('./ref_vCM.rds')

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "vCM" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("time_point", "celltype.lv4"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'time_point' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

#co-expression network analysis
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "v_CM", # the name of the group of interest in the group.by column
  group.by='celltype.lv4', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'SCT',
  layer = 'data'
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

pdf("Softpower_test.pdf")
# assemble with patchwork
wrap_plots(plot_list, ncol=2)

dev.off()

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'v_CM' # name of the topoligical overlap matrix written to disk
)
pdf("vCM_hdWGCNA_dendro.pdf")
PlotDendrogram(seurat_obj, main='v_CM hdWGCNA Dendrogram')
dev.off()

# need to run ScaleData first or else harmony throws an error:
#seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj
)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'celltype.lv4', group_name = 'v_CM'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "vCM"
)
pdf("vCM_kME.pdf")
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)
p
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)


# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)


# set new colors
new_colors <- list('vCM1' = '#073960', 'vCM2' = '#F94B4E', 'vCM3' ='#A84D4F', 'vCM4'='#1C80B2')
seurat_obj <- ResetModuleColors(seurat_obj, new_colors)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
pdf("vCM_module_featureplot_recolor.pdf")
# stitch together with patchwork
wrap_plots(plot_list, ncol=1)
dev.off()


ModuleNetworkPlot(
    seurat_obj, 
    outdir='ModuleNetworks_recolor', # new folder name
    n_inner = 10, # number of genes in inner ring
    n_outer = 20, # number of genes in outer ring
    n_conns = Inf, # show all of the connections
    plot_size=c(10,10), # larger plotting area
    vertex.label.cex=1 # font size
)

pdf("vCM_module_hub_network.pdf")
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
pdf("vCM_module_hub_network_umap.pdf")

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color, # color each point by WGCNA module
   size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

dev.off()

pdf("vCM_module_hub_network_umap_recolor.pdf")

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.1,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=3 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE

)
dev.off()


###Figure 5C
module_gene<-read.csv("vCM_module_assignment.csv", header = T)
module_gene.vCM1<- filter(module_gene, module_gene$module == "vCM1")
module_gene.vCM1 <- module_gene.vCM1[order(-module_gene.vCM1$kME_vCM1), ]
module_gene.vCM2<- filter(module_gene, module_gene$module == "vCM2")
module_gene.vCM2 <- module_gene.vCM2[order(-module_gene.vCM2$kME_vCM2), ]
module_gene.vCM3<- filter(module_gene, module_gene$module == "vCM3")
module_gene.vCM3 <- module_gene.vCM3[order(-module_gene.vCM3$kME_vCM3), ]
module_gene.vCM4<- filter(module_gene, module_gene$module == "vCM4")
module_gene.vCM4 <- module_gene.vCM4[order(-module_gene.vCM4$kME_vCM4), ]
module.gene<-list()
module.gene<-list("vCM1" = as.vector(module_gene.vCM1$gene_name[1:100]), 
                  "vCM2" = as.vector(module_gene.vCM2$gene_name[1:100]),
                  "vCM3" = as.vector(module_gene.vCM3$gene_name[1:100]),
                  "vCM4" = as.vector(module_gene.vCM4$gene_name[1:100]))
result <- compareCluster(geneCluster = module.gene, fun = enrichGO, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")

# Assuming you already have enrichGO result saved in `ego`
# Filter GO terms with padj < 0.05
go_result <- result@compareClusterResult %>%
  filter(p.adjust < 0.05)  # Filter GO terms with padj < 0.05
vCM1.go<-go_result[(go_result$Cluster == "vCM1"),]

# Add the -log10(p.adjust) to the dataframe
vCM1.go <- vCM1.go %>%
  mutate(log_padj = -log10(p.adjust))  # Add new column with -log10(p.adjust)

vCM1.go <- vCM1.go %>%
  mutate(
    GeneRatio_first = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)),  # Extract first number
    GeneRatio_second = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)), # Extract second number
    GeneRatio_value = GeneRatio_first / GeneRatio_second  # Calculate the actual ratio
  )


# Add labels for selected GO terms
select.go<-c("GO:0060047","GO:0042692","GO:0019932","GO:0009150","GO:0071805","GO:0006637")
select_go<-vCM1.go[(vCM1.go$ID %in% select.go),]
# Create a bubble plot
bubble_plot[["vCM1"]] <- ggplot(select_go, aes(x = GeneRatio_value, y = log_padj)) +
  geom_point(aes(size = Count), color = "#073960") +  # Bubble size is Count, color is padj
  scale_size(range = c(5, 10)) +  # Adjust bubble size range
  ylim(0, max(select_go$log_padj)) +
  labs(
    title = "GO Terms",
    x = "GeneRatio",
    y = "log_padj",
    size = "Gene Count",
    color = "Adjusted p-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1) )  # Rotate x-axis labels
bubble_plot[["vCM1"]]  <- bubble_plot[["vCM1"]]  + 
  geom_text_repel(
    data = select_go,  # Only label the top 5 GO terms
    aes(label = ID), 
    size = 5,  # Adjust text size for labels
    box.padding = 0.5, 
    point.padding = 0.5,
    segment.color = "grey50",  # Line color between the label and the dot
    segment.size = 0.5         # Line thickness
  )
# Print the bubble plot
pdf("../../../../Manuscript/Figure 4/vCM1_GO_result.pdf", height = 10, width = 12)
print(bubble_plot[["vCM1"]] )
dev.off()

vCM2.go<-go_result[(go_result$Cluster == "vCM2"),]

# Add the -log10(p.adjust) to the dataframe
vCM2.go <- vCM2.go %>%
  mutate(log_padj = -log10(p.adjust))  # Add new column with -log10(p.adjust)

vCM2.go <- vCM2.go %>%
  mutate(
    GeneRatio_first = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)),  # Extract first number
    GeneRatio_second = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)), # Extract second number
    GeneRatio_value = GeneRatio_first / GeneRatio_second  # Calculate the actual ratio
  )

# Add labels for selected GO terms
select.go<-c("GO:0106030","GO:0030111","GO:0048639","GO:0050769","GO:0003215","GO:0048738")
select_go<-vCM2.go[(vCM2.go$ID %in% select.go),]
# Create a bubble plot
bubble_plot[["vCM2"]]  <- ggplot(select_go, aes(x = GeneRatio_value, y = log_padj)) +
  geom_point(aes(size = Count), color = "#f94b4e") +  # Bubble size is Count, color is padj
  scale_size(range = c(5, 10)) +  # Adjust bubble size range
  ylim(0, max(select_go$log_padj)) +
  labs(
    title = "GO Terms",
    x = "GeneRatio",
    y = "log_padj",
    size = "Gene Count",
    color = "Adjusted p-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1) )  # Rotate x-axis labels
bubble_plot[["vCM2"]]  <- bubble_plot[["vCM2"]]  + 
  geom_text_repel(
    data = select_go,  # Only label the top 5 GO terms
    aes(label = ID), 
    size = 5,  # Adjust text size for labels
    box.padding = 0.5, 
    point.padding = 0.5,
    segment.color = "grey50",  # Line color between the label and the dot
    segment.size = 0.5         # Line thickness
  )

# Print the bubble plot
pdf("../../../../Manuscript/Figure 4/vCM2_GO_result.pdf", height = 10, width =8)
print(bubble_plot)
dev.off()

vCM4.go<-go_result[(go_result$Cluster == "vCM4"),]

# Add the -log10(p.adjust) to the dataframe
vCM4.go <- vCM4.go %>%
  mutate(log_padj = -log10(p.adjust))  # Add new column with -log10(p.adjust)

vCM4.go <- vCM4.go %>%
  mutate(
    GeneRatio_first = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)),  # Extract first number
    GeneRatio_second = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)), # Extract second number
    GeneRatio_value = GeneRatio_first / GeneRatio_second  # Calculate the actual ratio
  )

# Add labels for selected GO terms
select.go<-c("GO:0048738","GO:0055007","GO:0042391","GO:0001508","GO:0070374","GO:1902001","GO:0045214")
select_go<-vCM4.go[(vCM4.go$ID %in% select.go),]
# Create a bubble plot
bubble_plot[["vCM4"]]  <- ggplot(select_go, aes(x = GeneRatio_value, y = log_padj)) +
  geom_point(aes(size = Count), color = "#2b84b4") +  # Bubble size is Count, color is padj
  scale_size(range = c(5, 10)) +  # Adjust bubble size range
  ylim(0, max(select_go$log_padj)) +
  labs(
    title = "GO Terms",
    x = "GeneRatio",
    y = "log_padj",
    size = "Gene Count",
    color = "Adjusted p-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1) )  # Rotate x-axis labels
bubble_plot[["vCM4"]]  <- bubble_plot[["vCM4"]]  + 
  geom_text_repel(
    data = select_go,  # Only label the top 5 GO terms
    aes(label = ID), 
    size = 5,  # Adjust text size for labels
    box.padding = 0.5, 
    point.padding = 0.5,
    segment.color = "grey50",  # Line color between the label and the dot
    segment.size = 0.5         # Line thickness
  )

combined_plot <- plot_grid(bubble_plot[["vCM1"]], bubble_plot[["vCM2"]], bubble_plot[["vCM4"]], ncol = 3)

# Print the bubble plot
pdf("../../../../Manuscript/Figure 4/GO_result.pdf", height = 10, width = 15)
print(combined_plot)
dev.off()
