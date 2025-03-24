# ------------------------------------------------------------------------------
# Title: Preparation of Cytospace Input Files from Seurat Objects
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script generates input files for Cytospace from Seurat objects. 
# It processes single-cell and spatial transcriptomics (ST) data at different 
# time points (p0, p7, p14, p21). The script reads a reference Seurat object, 
# subsets it by time point, and exports formatted files required for Cytospace analysis.
#
# Dependencies:
# - Seurat
# - Cytospace Github Repo: https://github.com/digitalcytometry/cytospace
#
# Usage:
# Run this script in an R environment where the Seurat object is preloaded.
#
# ------------------------------------------------------------------------------


source('./cytospace/cytospace/Prepare_input_files/generate_cytospace_from_seurat_object.R')

ref<-readRDS("./ref_lv5.rds")
Idents(ref) <- 'celltype.lv2'
ref <- subset(ref, )
ref_p0 <- subset(ref, subset = time_point == "p0")
ref_p7 <- subset(ref, subset = time_point == "p7")
ref_p14 <- subset(ref, subset = time_point == "p14")
ref_p21 <- subset(ref, subset = time_point == "p21")

generate_cytospace_from_scRNA_seurat_object(ref, dir_out='./cytospace/', 
                                            fout_prefix='', write_sparse=FALSE, rna_assay='RNA')
generate_cytospace_from_ST_seurat_object(xenium.p0, dir_out='./cytospace/', 
                                         fout_prefix='p0_', write_sparse=FALSE, slice='slice1')
generate_cytospace_from_ST_seurat_object(xenium.p7, dir_out='./cytospace/', 
                                         fout_prefix='p7_', write_sparse=FALSE, slice='slice1')
generate_cytospace_from_ST_seurat_object(xenium.p14, dir_out='./cytospace/', 
                                         fout_prefix='p14_', write_sparse=FALSE, slice='slice1')
generate_cytospace_from_ST_seurat_object(xenium.p21, dir_out='./cytospace/', 
                                         fout_prefix='p21_', write_sparse=FALSE, slice='slice1')


generate_cytospace_from_scRNA_seurat_object(ref_p0, dir_out='./cytospace/lv2/full', 
                                            fout_prefix='p0_', write_sparse=FALSE, rna_assay='RNA')
generate_cytospace_from_scRNA_seurat_object(ref_p7, dir_out='./cytospace/lv2/full', 
                                            fout_prefix='p7_', write_sparse=FALSE, rna_assay='RNA')
generate_cytospace_from_scRNA_seurat_object(ref_p14, dir_out='./cytospace/lv2/full', 
                                            fout_prefix='p14_', write_sparse=FALSE, rna_assay='RNA')
generate_cytospace_from_scRNA_seurat_object(ref_p21, dir_out='./cytospace/lv2/full', 
                                            fout_prefix='p21_', write_sparse=FALSE, rna_assay='RNA')

