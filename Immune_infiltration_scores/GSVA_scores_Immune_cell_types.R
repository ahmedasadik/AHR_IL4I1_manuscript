#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript/AHR_signature/
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the generation of the immune infiltration GSVA scores for the different cell types
## across all TCGA tumors
################################################

## Source the functions and parameters files
source("./functions_and_parameters.R")

## Load libraries
library(GSVA)

## TCGA_voom
TCGA_voom <- readRDS("./Zenodo_download/TCGA_DGE_voom_annot.rds")

# Read the infiltration gene lists
cell_types_metagene_lists <- readRDS("./Zenodo_download/IPS_cell_types_gene_lists.rds")

## estimate the GSVA scores for the different cell populations
GSVA_voom_Immune_cells <- map(TCGA_voom, function(a,b,c){
  glists <- map(b, function(x)x$metagene)
  gsva(expr = t(a$E), glists, parallel.sz=c, method="gsva")
}, b=cell_types_metagene_lists, c=no_cores)

saveRDS(GSVA_voom_Immune_cells,"./Results/RDS/TCGA_GSVA_voom_28_cell_types_IPS.rds")
