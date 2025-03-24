# ------------------------------------------------------------------------------
# Title: Post-Cytospace Analysis and Visualization of Cell-Type Assignments
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script processes the output from Cytospace, and performs correlation 
# analysis between measured (Xenium) and imputed (Cytospace) expression values.
#
# Key Functions:
# - Loads assigned cell-type annotations from Cytospace outputs
# - Merges imputed gene expression data into the Xenium Seurat object
# - Computes Pearson and Spearman correlation distributions
# - Generates density plots of correlation values across different time points
#
# Dependencies:
# - Seurat (v5 compatibility enabled)
# - ggplot2, dplyr
#
# Output:
# - Correlation results saved as CSV files
# - Density plots exported as a PDF (Correlations.pdf)
#
# Usage:
# - Run in an R environment where the Xenium Seurat objects and Cytospace 
#   output files are accessible.
# ------------------------------------------------------------------------------



######### New Cytospace with Seurat Celltype assigned #######
setwd("./cytoSPACE/0824_run/output")

xenium.p0 <- readRDS("./robjs/xenium.p0.final.rds")
xenium.p7 <- readRDS("./robjs/xenium.p7.final.rds")
xenium.p14 <- readRDS("./robjs/xenium.p14.final.rds")
xenium.p21 <- readRDS("./robjs/xenium.p21.final.rds")

assigned_locations <- read.csv("singlecell_withoutSTcelltype_p0/assigned_locations.csv")
expression_data_path <- "./singlecell_withoutSTcelltype_p0/assigned_expression" 
options(Seurat.object.assay.version = "v5")

assigned_expression <- Read10X(data.dir = expression_data_path)

matching_rows <- assigned_locations[assigned_locations$SpotID %in% colnames(xenium.p0), ]
matching_rows <- matching_rows[match(colnames(xenium.p0), matching_rows$SpotID), ]

rownames(matching_rows) <- matching_rows$SpotID

xenium.p0$cytoSPACE.final.anno <- matching_rows[,"CellType"]
xenium.p0$nCount_cytoSPACE <- NULL
xenium.p0$nFeature_cytoSPACE <- NULL
xenium.p0$OriginalCID_snRNA <- matching_rows[,"OriginalCID"]
print(dim(assigned_expression))  # Print the dimensions of the matrix
print(head(rownames(assigned_expression)))  # Print a few gene names
print(head(colnames(assigned_expression)))# Print a few barcodes (UniqueCID)
head(assigned_locations$UniqueCID)

matching_rows <- assigned_locations[match(colnames(xenium.p0), assigned_locations$SpotID), ]
ordered_expression <- assigned_expression[, match(matching_rows$UniqueCID, colnames(assigned_expression))]
xenium.p0[["cytoSPACE"]] <- CreateAssay5Object(counts = ordered_expression)
saveRDS(xenium.p0, "./robjs/xenium.p0.0825.rds")


assigned_locations <- read.csv("singlecell_withoutSTcelltype_p7/assigned_locations.csv")
expression_data_path <- "./singlecell_withoutSTcelltype_p7/assigned_expression" 
options(Seurat.object.assay.version = "v5")

assigned_expression <- Read10X(data.dir = expression_data_path)
matching_rows <- assigned_locations[assigned_locations$SpotID %in% colnames(xenium.p7), ]
matching_rows <- matching_rows[match(colnames(xenium.p7), matching_rows$SpotID), ]
rownames(matching_rows) <- matching_rows$SpotID
xenium.p7$cytoSPACE.final.anno <- matching_rows[,"CellType"]
xenium.p7$nCount_cytoSPACE <- NULL
xenium.p7$nFeature_cytoSPACE <- NULL
xenium.p7$OriginalCID_snRNA <- matching_rows[,"OriginalCID"]

matching_rows <- assigned_locations[match(colnames(xenium.p7), assigned_locations$SpotID), ]
ordered_expression <- assigned_expression[, match(matching_rows$UniqueCID, colnames(assigned_expression))]
xenium.p7[["cytoSPACE"]] <- CreateAssay5Object(counts = ordered_expression)
saveRDS(xenium.p7, "./robjs/xenium.p7.0825.rds")




setwd("./cytoSPACE/0824_run/output")
assigned_locations <- read.csv("singlecell_withoutSTcelltype_p14/assigned_locations.csv")
expression_data_path <- "./singlecell_withoutSTcelltype_p14/assigned_expression" 
options(Seurat.object.assay.version = "v5")

assigned_expression <- Read10X(data.dir = expression_data_path)
matching_rows <- assigned_locations[assigned_locations$SpotID %in% colnames(xenium.p14), ]
matching_rows <- matching_rows[match(colnames(xenium.p14), matching_rows$SpotID), ]
rownames(matching_rows) <- matching_rows$SpotID
xenium.p14$cytoSPACE.final.anno <- matching_rows[,"CellType"]
xenium.p14$nCount_cytoSPACE <- NULL
xenium.p14$nFeature_cytoSPACE <- NULL
xenium.p14$OriginalCID_snRNA <- matching_rows[,"OriginalCID"]

matching_rows <- assigned_locations[match(colnames(xenium.p14), assigned_locations$SpotID), ]
ordered_expression <- assigned_expression[, match(matching_rows$UniqueCID, colnames(assigned_expression))]
xenium.p14[["cytoSPACE"]] <- CreateAssay5Object(counts = ordered_expression)
saveRDS(xenium.p14, "./robjs/xenium.p14.0825.rds")


setwd("./cytoSPACE/output")
assigned_locations <- read.csv("singlecell_withoutSTcelltype_p21/assigned_locations.csv")
expression_data_path <- "./singlecell_withoutSTcelltype_p21/assigned_expression" 
options(Seurat.object.assay.version = "v5")

assigned_expression <- Read10X(data.dir = expression_data_path)
matching_rows <- assigned_locations[assigned_locations$SpotID %in% colnames(xenium.p21), ]
matching_rows <- matching_rows[match(colnames(xenium.p21), matching_rows$SpotID), ]
rownames(matching_rows) <- matching_rows$SpotID
xenium.p21$cytoSPACE.final.anno <- matching_rows[,"CellType"]
xenium.p21$nCount_cytoSPACE <- NULL
xenium.p21$nFeature_cytoSPACE <- NULL
xenium.p21$OriginalCID_snRNA <- matching_rows[,"OriginalCID"]

matching_rows <- assigned_locations[match(colnames(xenium.p21), assigned_locations$SpotID), ]
ordered_expression <- assigned_expression[, match(matching_rows$UniqueCID, colnames(assigned_expression))]
xenium.p21[["cytoSPACE"]] <- CreateAssay5Object(counts = ordered_expression)
saveRDS(xenium.p21, "./robjs/xenium.p21.0825.rds")


##### Person correlation ##### 
library(ggplot2)
library(dplyr)

# Define a list of Seurat objects and their corresponding timepoint names
timepoint_objects <- list(xenium.p0 = xenium.p0, xenium.p7 = xenium.p7, 
                          xenium.p14 = xenium.p14, xenium.p21 = xenium.p21)
timepoint_names <- names(timepoint_objects)

# Function to calculate correlations and add a timepoint label
plot_correlation_distribution <- function(timepoint, timepoint_name, method = "pearson") {
  
  # Step 2: Identify the common genes between cytoSPACE and Xenium
  common_genes <- intersect(rownames(timepoint[["Xenium"]]), rownames(timepoint[["cytoSPACE"]]))
  
  # Step 3: Extract the measured (Xenium) and imputed (cytoSPACE) data for common genes
  measured_data <- subset(timepoint[["Xenium"]], features = common_genes)
  imputed_data <- subset(timepoint[["cytoSPACE"]], features = common_genes)
  
  # Normalize and scale the data
  measured_data <- NormalizeData(measured_data)
  # measured_data <- ScaleData(measured_data)
  imputed_data <- NormalizeData(imputed_data)
  # imputed_data <- ScaleData(imputed_data)
  
  # Step 4: Extract the assay data for the common genes before the loop to speed up calculations
  measured_data_matrix <- GetAssayData(measured_data, slot = "data")[common_genes, ]
  imputed_data_matrix <- GetAssayData(imputed_data, slot = "data")[common_genes, ]
  
  # Step 5: Randomly select 2,000 cells from the dataset
  total_cells <- ncol(measured_data_matrix)
  selected_cells <- sample(1:total_cells, 2000, replace = FALSE)
  
  # Step 6: Initialize a vector to store the correlation for the selected cells
  correlations <- numeric(length(selected_cells))
  
  # Step 7: Loop through each randomly selected cell and calculate the correlation based on the chosen method
  for (i in seq_along(selected_cells)) {
    cell_idx <- selected_cells[i]
    correlations[i] <- cor(measured_data_matrix[, cell_idx], 
                           imputed_data_matrix[, cell_idx], 
                           method = method)
  }
  
  # Step 8: Convert correlations to a data frame and add a timepoint label
  correlation_df <- data.frame(Correlation = correlations, Timepoint = timepoint_name)
  
  # Step 9: Return the correlation data frame
  return(correlation_df)
}

# Collect correlation data for all four time points
all_timepoints_df <- bind_rows(lapply(timepoint_names, function(timepoint_name) {
  plot_correlation_distribution(timepoint_objects[[timepoint_name]], timepoint_name, method = "pearson")
}))

# Plot density plots for all time points in one figure
ggplot(all_timepoints_df, aes(x = Correlation, fill = Timepoint)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Pearson Correlations Across Timepoints (Random 2000 Cells)",
       x = "Pearson Correlation", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("xenium.p0" = "blue", "xenium.p7" = "green", 
                               "xenium.p14" = "red", "xenium.p21" = "purple"))



all_timepoints_df_spearman <- bind_rows(lapply(timepoint_names, function(timepoint_name) {
  plot_correlation_distribution(timepoint_objects[[timepoint_name]], timepoint_name, method = "spearman")
}))
ggplot(all_timepoints_df_spearman, aes(x = Correlation, fill = Timepoint)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Spearman Correlations Across Timepoints (Random 2000 Cells)",
       x = "Spearman Correlation", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("xenium.p0" = "blue", "xenium.p7" = "green", 
                               "xenium.p14" = "red", "xenium.p21" = "purple"))



#### extract the genes name too ###
plot_correlation_distribution <- function(timepoint, timepoint_name, method = "pearson") {
  
  # Step 2: Identify the common genes between cytoSPACE and Xenium
  common_genes <- intersect(rownames(timepoint[["Xenium"]]), rownames(timepoint[["cytoSPACE"]]))
  
  # Step 3: Extract the measured (Xenium) and imputed (cytoSPACE) data for common genes
  measured_data <- subset(timepoint[["Xenium"]], features = common_genes)
  imputed_data <- subset(timepoint[["cytoSPACE"]], features = common_genes)
  
  # Normalize and scale the data
  measured_data <- NormalizeData(measured_data)
  #measured_data <- ScaleData(measured_data)
  imputed_data <- NormalizeData(imputed_data)
  #imputed_data <- ScaleData(imputed_data)
  
  # Step 4: Extract the assay data for the common genes
  measured_data_matrix <- GetAssayData(measured_data, slot = "data")[common_genes, ]
  imputed_data_matrix <- GetAssayData(imputed_data, slot = "data")[common_genes, ]
  
  # Step 5: Initialize a vector to store the correlation for each gene
  correlations <- numeric(length(common_genes))
  
  # Step 6: Calculate the correlation for each gene across all cells
  for (i in seq_along(common_genes)) {
    correlations[i] <- cor(measured_data_matrix[i, ], imputed_data_matrix[i, ], method = method)
  }
  
  # Step 7: Create a data frame for gene correlations
  correlation_df <- data.frame(Gene = common_genes, Correlation = correlations, Timepoint = timepoint_name)
  
  # Step 8: Find the top 10 genes with the highest correlation
  #top_genes_df <- correlation_df %>% arrange(desc(Correlation)) %>% head(1000)
  
  # Step 9: Return the data frame with the top genes for this timepoint
  return(correlation_df)
}
timepoint_objects <- list(xenium.p0 = xenium.p0, xenium.p7 = xenium.p7, 
                          xenium.p14 = xenium.p14, xenium.p21 = xenium.p21)
timepoint_names <- names(timepoint_objects)

# spearman
all_timepoints_df_spearman <- bind_rows(lapply(timepoint_names, function(timepoint_name) {
  plot_correlation_distribution(timepoint_objects[[timepoint_name]], timepoint_name, method = "spearman")
}))
write.csv(all_timepoints_df_spearman, "spearman_correlation_across_timepoints.csv", row.names = FALSE)
all_timepoints_df_spearman<- read.csv("./spearman_correlation_across_timepoints.csv")

ggplot(all_timepoints_df_spearman, aes(x = Correlation, fill = Timepoint)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Spearman Correlations Across Timepoints",
       x = "Spearman Correlation", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("xenium.p0" = "blue", "xenium.p7" = "green", 
                               "xenium.p14" = "red", "xenium.p21" = "purple"))

all_timepoints_df <- all_timepoints_df %>% arrange(desc(Correlation)) 


# pearson
all_timepoints_df <- bind_rows(lapply(timepoint_names, function(timepoint_name) {
  plot_correlation_distribution(timepoint_objects[[timepoint_name]], timepoint_name, method = "pearson")
}))
write.csv(all_timepoints_df, "pearson_correlation_across_timepoints.csv", row.names = FALSE)

all_timepoints_df<- read.csv("pearson_correlation_across_timepoints.csv")



p1<- ggplot(all_timepoints_df, aes(x = Correlation, fill = Timepoint)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Pearson Correlations Across Timepoints",
       x = "Pearson Correlation", y = "Density") +
  theme_classic() +
  scale_fill_manual(values = c("xenium.p0" = "#FDF6B5" , "xenium.p7" = "#F9C480", 
                               "xenium.p14" = "#F28D6F", "xenium.p21" ="#E24C80"))

p2 <- ggplot(all_timepoints_df_spearman, aes(x = Correlation, fill = Timepoint)) +
  geom_density(alpha = 0.5) +
  labs(title = "Comparison of Spearman Correlations Across Timepoints",
       x = "Spearman Correlation", y = "Density") +
  theme_classic() +
  scale_fill_manual(values = c("xenium.p0" = "#FDF6B5" , "xenium.p7" = "#F9C480", 
                               "xenium.p14" = "#F28D6F", "xenium.p21" ="#E24C80"))
pdf("../plot/Correlations.pdf")
print(p1)
print(p2)
dev.off()



