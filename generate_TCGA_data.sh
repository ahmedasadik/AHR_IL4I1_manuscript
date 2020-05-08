#!/bin/bash

# A script to setup the raw data objects
# These do not need to be run if you download processed data from Zenodo

# 1: Generating_TCGA_data_objects
cd Generating_TCGA_data_objects
./Generating_TCGA_WGCNA_objects.R
./Generating_TCGA_data_objects.R

# 2: Generating_TCGA_immune_infiltration_scores
cd Generating_TCGA_immune_infiltration_scores
./GSVA_scores_Immune_cell_types.R
./Infilteration_gene_signatures.R
