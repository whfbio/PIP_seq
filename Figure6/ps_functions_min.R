# ============================================================
# Helper functions for Perturb-seq module-score analysis
# (minimal subset extracted from original functions.R)
# ============================================================

#' Count cells expressing a gene above a threshold
#' @param object Seurat object
#' @param genes Character, gene/feature name
#' @param threshold Numeric, expression threshold on raw counts
#' @return Integer count of cells with counts > threshold
get_expressing_cells_num <- function(object, genes, threshold = 3) {
  counts <- object[["RNA"]]@layers$counts
  colnames(counts) <- colnames(object)
  rownames(counts) <- rownames(object)

  if (!genes %in% rownames(counts)) return(NA_integer_)
  sum(counts[genes, ] > threshold)
}

#' Cells expressing both features above a threshold
get_double_expressing_cells <- function(object, gene1, gene2, threshold = 2) {
  counts <- object[["RNA"]]@layers$counts
  colnames(counts) <- colnames(object)
  rownames(counts) <- rownames(object)

  if (!(gene1 %in% rownames(counts) && gene2 %in% rownames(counts))) return(NA)
  keep <- (counts[gene1, ] > threshold) & (counts[gene2, ] > threshold)
  colnames(counts)[keep]
}

get_double_expressing_cells_num <- function(object, gene1, gene2, threshold = 2) {
  cells <- get_double_expressing_cells(object, gene1, gene2, threshold = threshold)
  if (is.na(cells)[1]) return(NA_integer_)
  length(cells)
}

#' Cells expressing either feature above a threshold (used for QC summaries)
get_single_expressing_cells_num <- function(object, gene1, gene2, threshold = 2) {
  counts <- object[["RNA"]]@layers$counts
  colnames(counts) <- colnames(object)
  rownames(counts) <- rownames(object)

  if (!(gene1 %in% rownames(counts) && gene2 %in% rownames(counts))) return(NA_integer_)
  sum((counts[gene1, ] > threshold) | (counts[gene2, ] > threshold))
}

#' Join TRUE guide labels for a cell (row of logicals)
get_true_sgRNA_names <- function(row) {
  sgRNA_names <- names(row[row])
  paste(sgRNA_names, collapse = ", ")
}

#' Compute a custom module score from scaled expression (Seurat v5 layers).
#' This reproduces the behavior in the original script:
#'  - z-score each gene across cells (on the scaled layer), then average across genes
#'  - cap scores to quantile cutoffs to reduce outlier influence
#'
#' @param genes Character vector of genes in the module
#' @param seurat Seurat object
#' @param cutoff Numeric length-2, lower/upper quantiles used to cap scores
#' @param assay Assay name
#' @return Numeric vector of scores (length = number of cells)
module_score <- function(genes, seurat, cutoff = c(.005, .995), assay = "RNA") {
  genes <- unlist(genes)
  score <- numeric(ncol(seurat))
  gene_num <- 0

  # Seurat v5: scaled data is stored in layers$scale.data
  expression_data <- as.data.frame(seurat[[assay]]@layers$scale.data)
  rownames(expression_data) <- rownames(seurat@assays[[assay]])
  colnames(expression_data) <- colnames(seurat@assays[[assay]])

  for (gene in genes) {
    if (!gene %in% rownames(expression_data)) next
    gene_data <- data.frame(expression_data[gene, ])
    if (all(gene_data == 0, na.rm = TRUE)) next

    gene_num <- gene_num + 1
    row_means <- rowMeans(gene_data)
    row_sd <- apply(gene_data, 1, sd)

    standardized <- t(apply(gene_data, 1, function(row) (row - row_means) / row_sd))
    score <- score + standardized

    if (any(is.nan(score))) stop("NaN encountered in module score calculation.")
  }

  if (gene_num == 0) return(rep(0, length(score)))

  score <- score / gene_num
  thresh <- quantile(score, cutoff)
  score[score < thresh[1]] <- thresh[1]
  score[score > thresh[2]] <- thresh[2]
  score
}
