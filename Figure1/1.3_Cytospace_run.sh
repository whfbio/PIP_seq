#!/bin/bash
# ------------------------------------------------------------------------------
# Title: Running Cytospace Alignment on Spatial Transcriptomics and scRNA-seq Data
# Author: Yiran Song
# Date: March 18, 2025
# Description:
# This script automates the execution of Cytospace using single-cell RNA-seq 
# and spatial transcriptomics data across multiple time points (p0, p7, p14, p21).
# It runs Cytospace on a high-performance computing (HPC) cluster using sbatch.
#
# Dependencies:
# - Cytospace (requires conda environment)
#
# Usage:
# Run this script in a terminal with the appropriate SLURM queue setup.
# ------------------------------------------------------------------------------


declare -a time_points=("p0" "p7" "p14" "p21")

module add anaconda/2021.11
conda activate cytospace_v1.1.0

for time_point in "${time_points[@]}"
do
    output_folder="output/singlecell_withoutSTcelltype_${time_point}"
    mkdir -p $output_folder
    stc_path="${time_point}_st_cell_type_labels.txt"
    st_path="${time_point}_ST_data.txt"
    scRNA_path="${time_point}_scRNA_data.txt" 
    cell_type_path="${time_point}_cell_type_labels.txt" 
    coordinates_path="${time_point}_Coordinates.txt"
    job_name="cyto_${time_point}"

    sbatch --job-name=$job_name \
           --output=${job_name}_%j.txt \
           --time=3-00:00:00 \
           --mem=700G \
           --cpus-per-task=16 \
           --wrap="cytospace --single-cell --scRNA-path $scRNA_path --cell-type-path $cell_type_path --st-path $st_path --st-cell-type-path $stc_path --coordinates-path $coordinates_path --output-folder $output_folder --number-of-processors 16 --number-of-selected-spots 10000"
done