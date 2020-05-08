#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script show how IDO1 and TDO2 strongly associate with AHR activation
#################################################

## Load libraries
library(purrr)
library(edgeR)
library(corrplot)

## Source the functions and parameters files
source("../functions_and_parameters.R")

## Load Data
# TCGA_counts
TCGA_counts <- readRDS("../Zenodo_download/TCGA_counts.rds")

# TCGA_voom
TCGA_voom <- readRDS("../Zenodo_download/TCGA_DGE_voom_annot.rds")

# Read the AHR signature file
overlapping_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", sep = "\t", stringsAsFactors = F)

####################################################################################
## AHR activation is associated with TDO2 and IDO1 expression in different tumors ##
####################################################################################
# Use the TCGA voom data set to avoid additional calculation steps, then run gene set testing

# TDO2 groups separated by the median into high and low. The comparison is high-low.
TDO2_GSA_roast_00_incr <- map2(TCGA_counts, TCGA_voom, GSA_roast_hi_lo_grps_FUN, "TDO2", overlapping_genes$Gene, 0) %>% do.call(rbind,.)
write.table(TDO2_GSA_roast_00_incr,"../Results/Tables/TDO2_GSA_roast_00_incr.txt", sep = "\t")

# IDO1 groups separated by the median into high and low. The comparison is high-low.
IDO1_GSA_roast_00_incr <- map2(TCGA_counts, TCGA_voom, GSA_roast_hi_lo_grps_FUN, "IDO1", overlapping_genes$Gene, 0) %>% do.call(rbind,.)
write.table(IDO1_GSA_roast_00_incr,"../Results/Tables/IDO1_GSA_roast_00_incr.txt", sep = "\t")

## Represent the up and down regulation as percentages of the proportional genes that are significantly different
GSA_corrplot_r <- data.frame(IDO1_dn=IDO1_GSA_roast_00_incr$PropDown*-1,
                             IDO1_up=IDO1_GSA_roast_00_incr$PropUp,
                             blank=rep(0,32),
                             TDO2_dn=TDO2_GSA_roast_00_incr$PropDown*-1,
                             TDO2_up=TDO2_GSA_roast_00_incr$PropUp)

GSA_corrplot_p <- data.frame(IDO1_dn=IDO1_GSA_roast_00_incr$FDR,
                             IDO1_up=IDO1_GSA_roast_00_incr$FDR,
                             blank=rep(1,32),
                             TDO2_dn=TDO2_GSA_roast_00_incr$FDR,
                             TDO2_up=TDO2_GSA_roast_00_incr$FDR)

rownames(GSA_corrplot_r) <- rownames(GSA_corrplot_p) <- gsub("TCGA_","",tcga_names)
colnames(GSA_corrplot_r) <- colnames(GSA_corrplot_p) <- c("IDO1_low", "IDO1_high"," ", "TDO2_low", "TDO2_high")

pdf("../Results/Figures/GSA_IDO_TDO_median.pdf", width = 8, height = 12, pointsize = 12)
corrplot(corr = as.matrix(GSA_corrplot_r), p.mat = as.matrix(GSA_corrplot_p), insig = "blank", mar=c(0,0,1,0),
         cl.length = 3, col = colorRampPalette(colors = c("blue", "white", "red"))(11), addgrid.col = NA,
         tl.cex = 0.7, is.corr = F, method = "pie", cl.pos = "b", cl.lim = c(-1,1), outline = "black",
         tl.col = "black")
dev.off()