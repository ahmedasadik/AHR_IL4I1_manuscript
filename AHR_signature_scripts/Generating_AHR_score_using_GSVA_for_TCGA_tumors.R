#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This file explains the generation of a GSVA score for the AHR signature
#################################################

## Load libraries
library(edgeR)
library(GSVA)

## Source the functions and parameters files
source("../functions_and_parameters.R")

# Loading TCGA DGE object raw counts before normalization as DGElists
TCGA_DGE <- readRDS("../Zenodo_download/TCGA_DGE.rds")

# Read the AHR signature file
overlapping_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", sep = "\t", stringsAsFactors = F)

# Extracting the AHR ENSG_IDs
AHR_ensg <- TCGA_DGE$TCGA_ACC$genes[TCGA_DGE$TCGA_ACC$genes$external_gene_name %in% overlapping_genes$Gene,]

# generate the GSVA enrichment scores of all TCGA tumors
TCGA_gsva <- map(TCGA_DGE, safely(function(dge){(gsva(cpm(dge$counts, prior.count = 1, log = T),list(AHR_ensg$ensembl_gene_id), parallel.sz=no_cores))}))
saveRDS(TCGA_gsva, "../Results/RDS/TCGA_GSVA_scores_safely.rds")
