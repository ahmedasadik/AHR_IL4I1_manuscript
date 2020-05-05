#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script explains how WGCNA was performed for all TCGA tumors
#################################################

## Load libraries
library(purrr)
library(WGCNA)

## Source the functions and parameters files
source("./functions_and_parameters.R")

## load TCGA_voom
TCGA_DGE_voom_annots <- readRDS("./Zenodo_download/TCGA_DGE_voom_annots.rds")

# Enabling WGCNA threading
enableWGCNAThreads(nThreads = no_cores)

# Estimating soft threshold power
# Selecting the powers for soft thresholding
powers <- c(c(1:10), seq(from = 12, to=40, by=2))

# Select soft threshold
TCGA_sft <- map(TCGA_DGE_voom_annots, function(dge, pw_vec){
  pickSoftThreshold((dge$E), dataIsExpr = TRUE, networkType = "signed hybrid",
                    blockSize = dim(dge$E)[[2]], powerVector = pw_vec, verbose = 5)
}, powers)

saveRDS(TCGA_sft, "./Results/RDS/TCGA_sft.rds")

## RUN WGCNA and generate coexpression modules
TCGA_TOMS <- map(1:length(TCGA_DGE_voom_annots), function(i, dge_voom, sft_power, nc){
  blockwiseModules(datExpr = (dge_voom[[i]]$E),corType = "bicor",
                   maxBlockSize = dim(dge_voom[[i]]$E)[[2]], nThreads = nc,
                   randomSeed = 0861, power = sft_power[[i]]$powerEstimate,
                   networkType = "signed hybrid",
                   TOMType = "signed", TOMDenom = "mean", saveTOMs = TRUE,
                   saveTOMFileBase = names(dge_voom)[i],
                   pamStage = TRUE, pamRespectsDendro = FALSE,verbose = 3)},
  dge_voom=TCGA_DGE_voom_annots, sft_power=TCGA_sft, nc=no_cores)

saveRDS(TCGA_TOMS, "./Results/RDS/TCGA_TOMS_bicor.rds")

## Extract the gene names for the different modules
TCGA_modules_genes <- map2(TCGA_TOMS, TCGA_DGE_voom_annots, function(tom, dge){
  df <- data.frame(genes=colnames(dge$E), module=tom$colors, stringsAsFactors = F)
})
saveRDS(TCGA_modules_genes, "./Results/RDS/TCGA_modules_genes.rds")

