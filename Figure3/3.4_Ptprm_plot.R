# ------------------------------------------------------------------------------
# Title: Plotting TF Weights for Target Gene (Meis1)
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script loads gene regulatory network (GRN) output from Commot 
# and visualizes the weight distribution of transcription factors (TFs)
# targeting "Meis1" across different time points (P0, P7, P14, P21).
#
# Key Functions:
# - Reads GRN output files
# - Filters TF interactions for target gene "Meis1"
# - Identifies top 5 TFs by weight
# - Generates scatter plots with labeled top TFs
#
# Dependencies:
# - ggplot2, ggrepel, readr, dplyr
#
# Output:
# - PDF plots showing TF-Target interactions across different time points
#
# Usage:
# - Run after GRN inference using Commot.
# ------------------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)

# Define time points
time_points <- c("p0", "p7", "p14", "p21")

# Function to process data and generate plots
plot_tf_weights <- function(time_point) {
  
  # Load gene regulatory network data
  file_path <- paste0("./grn/commot_grn_output_", time_point, "_cytospace_data.tsv")
  linkList <- read_delim(file_path, delim = "\t", col_names = FALSE)
  colnames(linkList) <- c('TF', 'Target', 'Weight')
  
  # Convert to appropriate data types
  linkList$TF <- as.character(linkList$TF)
  linkList$Target <- as.character(linkList$Target)
  linkList$Weight <- as.numeric(linkList$Weight)
  
  # Filter for target gene "Meis1"
  filtered_data <- linkList %>% filter(Target == "Meis1")
  
  # Identify top 5 TFs by weight
  top5 <- filtered_data %>% arrange(desc(Weight)) %>% head(5)
  filtered_data <- filtered_data %>%
    mutate(is_top5 = TF %in% top5$TF)
  
  # Generate plot
  p <- ggplot(filtered_data, aes(x = seq_along(Weight), y = Weight)) +
    geom_point(aes(color = is_top5), size = 3) +  
    geom_text_repel(aes(label = ifelse(is_top5, TF, "")),
                    size = 3, box.padding = 0.35, point.padding = 0.3,
                    segment.color = "grey50", max.overlaps = Inf) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +  
    labs(x = "Rank", y = "Weight", title = paste("TF Weights for Meis1 (", time_point, ")")) +
    theme_classic()
  
  # Save plot
  output_path <- paste0("ptprm_", time_point, ".pdf")
  ggsave(output_path, plot = p, width = 10, height = 6)
}

# Run the function for all time points
lapply(time_points, plot_tf_weights)
