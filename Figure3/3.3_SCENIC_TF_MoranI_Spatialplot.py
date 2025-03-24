# ------------------------------------------------------------------------------
# Title: SCENIC Transcription Factor Moran's I Analysis
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script calculates Moran's I scores for individual transcription factors (TFs)
# using spatial transcriptomics (Xenium) data and SCENIC-inferred activity.
# 
# Key Functions:
# - Loads SCENIC regulon activity into the AnnData object
# - Computes Moran's I for a predefined set of TFs
# - Generates spatial autocorrelation plots for each TF
#
# Dependencies:
# - anndata, scanpy, squidpy, numpy, pandas, matplotlib, seaborn
#
# Output:
# - Updated AnnData object with SCENIC activity (`xenium_total_SCENIC.h5ad`)
# - CSV file with Moran's I scores for all TFs
# - Plots of Moran's I scores saved as PDFs
#
# Usage:
# - Ensure SCENIC results are available (from `3.1_SCENIC_Integration.R`)
# ------------------------------------------------------------------------------


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc
import squidpy as sq
import os


adata = ad.read_h5ad("./xenium_total.h5ad")
adata.obs['combined_cell_id'] = adata.obs['time_point'].astype(str) + "_" + adata.obs['cell_id'].astype(str)
new_matrix = pd.read_feather("./scenic_mat_mapped_data.feather")
cell_ids = adata.obs['combined_cell_id'].values
matching_columns = new_matrix.columns.isin(cell_ids)
new_matrix_subset = new_matrix.loc[:, matching_columns]
print(f"Subset matrix shape: {new_matrix_subset.shape}")
new_matrix_subset_reordered = new_matrix_subset[adata.obs['combined_cell_id']]
new_matrix_subset_reordered_T = new_matrix_subset_reordered.T  # Shape (403168, 220)
print(new_matrix_subset_reordered_T.shape)
# Add the transposed matrix as a new entry in .obsm
adata.obsm['SCENIC'] = new_matrix_subset_reordered_T.values
adata.write("./xenium_total_SCENIC.h5ad")

# Caculate Moran'I score of TF
# List of transcription factors (TFs)
tf_names = [
    "Meis2", "Pbx1", "Irx5", "Irx2", "Sox4", "Rorb", "Bcl11a", "Gata6", "Mef2a", "Gata4",
    "Rb1", "Suz12", "Tfdp1", "Mybl2", "Brca1", "E2f1", "E2f2", "Mybl1", "E2f8", "Sf1",
    "Chd1", "E2f7", "E2f3", "Pbx3", "Hmga2", "Srebf2", "Sp1", "Gabpb1", "Crem", "Gtf3c2",
    "Smarca4", "Phf8", "Taf1", "Sin3a", "Arnt", "Rcor1", "Sp3", "Creb1", "Stat5a", "Ppard",
    "Bach1", "Nfe2l2", "Max", "Foxo3", "Tef", "Cebpb", "Dbp", "Dnmt3a", "Nr4a1", "Atf6",
    "Atf1", "Ppargc1a", "Hif1a", "Esrrg", "Rarb", "Mef2d", "Esrra", "Thrb", "Esrrb", "Ppara",
    "Stat3", "Atf3", "Fosl2", "Pax5", "Zfp831", "Batf", "Hsf3", "Mafb", "Runx3", "Irf4",
    "Irf5", "Arid3a", "Maf", "Jund", "Jun", "Tead2", "Setdb1", "E2f6", "Gabpa", "Mafg",
    "Zfp410", "Zfp397", "Rfx2", "Hlf", "Vdr", "Esr2", "Nfic", "Foxp4", "Hivep2", "Egr3",
    "Hivep1", "Pax3", "Sox10", "Gata3", "Smad9", "Lef1", "Gata5", "Nr2f2", "Rela", "Stat2",
    "Irf9", "Irf1", "Stat1", "Sox7", "Tal1", "Sox17", "Gata2", "Elk3", "Ets2", "Tbx3",
    "Mef2c", "Klf10", "Egr1", "Irf8", "Creb3l2", "Tcf7l2", "Tcf7l1", "Foxp2", "Gli3", "Gli2",
    "Gli1", "Nfkb1", "Zfp385a", "Yy1", "Kdm5a", "Irf2", "Nr3c1", "Klf2", "Klf4", "Klf3",
    "Zmiz1", "Foxo1", "Foxp1", "Etv6", "Fli1", "Erg", "Ets1", "Elk4", "Klf6", "Pou2f2",
    "Ikzf1", "Runx1", "Fez1", "Stat6", "Relb", "Cebpg", "Meis1-extended", "Tbx5-extended",
    "Xrcc4-extended", "Ezh2-extended", "Rad21-extended", "Tcf12-extended", "Zfp369-extended",
    "Brf1-extended", "Tbl1xr1-extended", "Nr2c2-extended", "Atf2-extended", "Vezf1-extended",
    "Zfp612-extended", "Cpeb1-extended", "Nr1d2-extended", "Klf15-extended", "Bcl6-extended",
    "Ddit3-extended", "Bhlhe40-extended", "Bdp1-extended", "Stat5b-extended", "Elf2-extended",
    "Cux1-extended", "Rreb1-extended", "Nr3c2-extended", "Rxrg-extended", "Cebpa-extended",
    "Lyl1-extended", "Hdac2-extended", "Ubtf-extended", "Zfp384-extended", "Zfp143-extended",
    "Usf2-extended", "Foxk1-extended", "Rest-extended", "Zfp110-extended", "Hoxd4-extended",
    "Smarcc2-extended", "Nfe2l1-extended", "Rbbp5-extended", "Zfp523-extended", "Etv5-extended",
    "Ep300-extended", "Nr1d1-extended", "Trp73-extended", "Isx-extended", "Tbx2-extended",
    "Nfatc4-extended", "Nr1h3-extended", "Tead4-extended", "Prdm1-extended", "Pparg-extended",
    "Nr5a2-extended", "Foxc1-extended", "Creb5-extended", "Nfat5-extended", "Twist2-extended",
    "Tbx18-extended", "Ar-extended", "Creb3l1-extended", "Tcf21-extended", "Mta3-extended",
    "Msc-extended", "Pknox2-extended", "Cebpd-extended", "Glis2-extended", "Ctcf-extended",
    "Mxi1-extended", "Elf1-extended", "Bclaf1-extended", "Smad1-extended", "Sox13-extended",
    "Tcf4-extended", "Zbtb7a-extended"
]

tf_name_to_index = {tf: i for i, tf in enumerate(tf_names)}
index_to_tf_name = {v: k for k, v in tf_name_to_index.items()}


# Caculate the Moran'I for each time point
# p0
subset_adata = adata[adata.obs['timepoint'] == 'p0']
sq.gr.spatial_autocorr(
    subset_adata,
    connectivity_key='spatial_connectivities',  # Uses the precomputed spatial connectivities
    genes=None,  # Since we are using obsm, not genes
    mode='moran',  # Calculate Moran's I
    attr='obsm',  # SCENIC data is stored in obsm
    layer='SCENIC',  # Refers to the key in obsm
    n_perms=100,  # Set to None for normality assumption; you can increase this for permutation tests
    corr_method='fdr_bh',  # Correction for multiple testing (FDR)
    transformation=True,  # Normalize spatial connectivities
    n_jobs=4,  # Number of parallel jobs (adjust based on your system)
    show_progress_bar=True
)
subset_adata.uns['moranI']
subset_adata.uns['moranI'].index = subset_adata.uns['moranI'].index.map(index_to_tf_name)
subset_adata.uns["moranI"].to_csv("./TF_moranI_p0.csv")
# p7
subset_adata = adata[adata.obs['timepoint'] == 'p7']
sq.gr.spatial_autocorr(
    subset_adata,
    connectivity_key='spatial_connectivities',  # Uses the precomputed spatial connectivities
    genes=None,  # Since we are using obsm, not genes
    mode='moran',  # Calculate Moran's I
    attr='obsm',  # SCENIC data is stored in obsm
    layer='SCENIC',  # Refers to the key in obsm
    n_perms=100,  # Set to None for normality assumption; you can increase this for permutation tests
    corr_method='fdr_bh',  # Correction for multiple testing (FDR)
    transformation=True,  # Normalize spatial connectivities
    n_jobs=4,  # Number of parallel jobs (adjust based on your system)
    show_progress_bar=True
)
subset_adata.uns["moranI"]
subset_adata.uns['moranI'].index = subset_adata.uns['moranI'].index.map(index_to_tf_name)
subset_adata.uns["moranI"].to_csv("./TF_moranI_p7.csv")
# p14
subset_adata = adata[adata.obs['timepoint'] == 'p14']
sq.gr.spatial_autocorr(
    subset_adata,
    connectivity_key='spatial_connectivities',  # Uses the precomputed spatial connectivities
    genes=None,  # Since we are using obsm, not genes
    mode='moran',  # Calculate Moran's I
    attr='obsm',  # SCENIC data is stored in obsm
    layer='SCENIC',  # Refers to the key in obsm
    n_perms=100,  # Set to None for normality assumption; you can increase this for permutation tests
    corr_method='fdr_bh',  # Correction for multiple testing (FDR)
    transformation=True,  # Normalize spatial connectivities
    n_jobs=4,  # Number of parallel jobs (adjust based on your system)
    show_progress_bar=True
)

subset_adata.uns['moranI'].index = subset_adata.uns['moranI'].index.map(index_to_tf_name)
subset_adata.uns["moranI"].to_csv("./TF_moranI_p14.csv")
subset_adata.uns["moranI"]
# p21
subset_adata = adata[adata.obs['timepoint'] == 'p21']
sq.gr.spatial_autocorr(
    subset_adata,
    connectivity_key='spatial_connectivities',  # Uses the precomputed spatial connectivities
    genes=None,  # Since we are using obsm, not genes
    mode='moran',  # Calculate Moran's I
    attr='obsm',  # SCENIC data is stored in obsm
    layer='SCENIC',  # Refers to the key in obsm
    n_perms=100,  # Set to None for normality assumption; you can increase this for permutation tests
    corr_method='fdr_bh',  # Correction for multiple testing (FDR)
    transformation=True,  # Normalize spatial connectivities
    n_jobs=4,  # Number of parallel jobs (adjust based on your system)
    show_progress_bar=True
)

subset_adata.uns['moranI'].index = subset_adata.uns['moranI'].index.map(index_to_tf_name)
subset_adata.uns["moranI"].to_csv("./TF_moranI_p21.csv")
subset_adata.uns["moranI"]


# Create individual TF plots, Figure 3D
output_folder = "figure_v2"
os.makedirs(output_folder, exist_ok=True)

# Loop over all columns in SCENIC (TFs)
for i, tf_name in enumerate(scenic_tf_names):
    # Add the TF score to obs with a prefix to avoid conflicts
    adata.obs[f"SCENIC_{tf_name}"] = adata.obsm['SCENIC'][:, i]
    
    # Calculate global min and max values for this TF across all time points
    global_min = adata.obs[f"SCENIC_{tf_name}"].min()
    global_max = adata.obs[f"SCENIC_{tf_name}"].max()

    # Create a new figure with 4 subplots (one for each time point)
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))  # One row, 4 columns for the 4 time points
    
    # Subset and plot for each time point
    for j, time_point in enumerate(['p0', 'p7', 'p14', 'p21']):
        # Subset the data for the current time point
        subset_adata = adata[adata.obs['timepoint'] == time_point].copy()
        
        # Plot using sc.pl.embedding with consistent vmin and vmax
        sc.pl.embedding(
            subset_adata,
            basis='spatial_fov',  # Use 'spatial' embedding for the plot
            color=f"SCENIC_{tf_name}",  # Color the plot by the TF with the prefix
            size=20,  # Adjust point size
            frameon=False,  # Remove the frame
            ax=axes[j],  # Use the subplot axis for each time point
            vmin=global_min,  # Set consistent minimum color scale
            vmax=global_max,  # Set consistent maximum color scale
            show=False,  # Do not show the plot in the notebook
            title=f"{tf_name} - {time_point}",  # Title showing the TF and time point
        )
    
    # Save each TF plot as a separate PDF in the figures folder
    plt.tight_layout()  # Adjust layout to avoid overlap
    output_path = os.path.join(output_folder, f"{tf_name}_spatial_embedding.pdf")
    plt.savefig(output_path, format='pdf')  # Save the plot for this TF
    plt.close()  # Close the figure to free up memory
