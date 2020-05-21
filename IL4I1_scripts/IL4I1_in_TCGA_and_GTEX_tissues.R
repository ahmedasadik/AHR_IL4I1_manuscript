#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This files describes the median expression of IDO1, TDO2, IDO2 and IL4I1 in TCGA and GTEX tissues
################################################

## Load libraries
library(purrr)
library(ggplot2)
library(corrplot)
library(reshape)
library(gridExtra)
library(ggbeeswarm)
library(ggpubr)

## Read TCGA_TPMs
TCGA_TPMs <- readRDS("../Zenodo_download/TCGA_TPMs.rds")

## estimate the mean or median of IDO1, IDO2, and TDO2 in the tumors and normal tissue
TCGA_enz_means <- map_df(TCGA_TPMs[-14], function(dge){
  IDO1_idx <- which(rownames(dge)=="IDO1")
  IDO2_idx <- which(rownames(dge)=="IDO2")
  TDO2_idx <- which(rownames(dge)=="TDO2")
  IL4I1_idx <- which(rownames(dge)=="IL4I1")
  tpm_m_sel <- dge[c(IDO1_idx, IDO2_idx, TDO2_idx, IL4I1_idx),] %>% t()
  colnames(tpm_m_sel) <- c("IDO1", "IDO2", "TDO2", "IL4I1")
  tpm_m_sel <- log2(tpm_m_sel)
  tpm_v_sel_med <- apply(tpm_m_sel, 2, median)
  names(tpm_v_sel_med) <- colnames(tpm_m_sel)
  tpm_v_sel_med
}) %>% as.data.frame() %>% t()

colnames(TCGA_enz_means) <- c("IDO1", "IDO2", "TDO2", "IL4I1")
TCGA_enz_means <- as.data.frame(TCGA_enz_means)
TCGA_enz_means[which(is.infinite(TCGA_enz_means$IDO2)),2] <- 0

rownames(TCGA_enz_means) <- gsub("TCGA_","", rownames(TCGA_enz_means))
pdf("../Results/Figures/bubbleheatmap_plot_IDO1_TDO2_IDO2_IL4I1_TCGA.pdf",width = 8, height = 12, pointsize = 24)
corrplot(corr = t(as.matrix(TCGA_enz_means)), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")
dev.off()

#################
## GTEX TPMs
gtex_tpms <- data.table::fread("../Zenodo_download/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")
gtex_annot <- read.delim("../Zenodo_download/GTEX_annot", sep = "\t", stringsAsFactors = F)

GTEX_IDO1_idx <- which(gtex_tpms$Description=="IDO1")
GTEX_IDO2_idx <- which(gtex_tpms$Description=="IDO2")
GTEX_TDO2_idx <- which(gtex_tpms$Description=="TDO2")
GTEX_IL4I1_idx <- which(gtex_tpms$Description=="IL4I1")

gtex_tpms_enz <- gtex_tpms[c(GTEX_IDO1_idx, GTEX_IDO2_idx, GTEX_TDO2_idx, GTEX_IL4I1_idx),] %>% t()

colnames(gtex_tpms_enz) <- c("IDO1", "IDO2", "TDO2", "IL4I1")
gtex_tpms_enz <- gtex_tpms_enz[-c(1:2),]
gtex_tpms_enz_log2 <- apply(gtex_tpms_enz, 2, as.numeric) %>% log2()
rownames(gtex_tpms_enz_log2) <- rownames(gtex_tpms_enz)

new_gtex_annot <- gtex_annot[match(rownames(gtex_tpms_enz_log2), gtex_annot$SAMPID)%>%.[!is.na(.)],]

gtex_counts <- as.data.frame(gtex_tpms_enz_log2)[match(new_gtex_annot$SAMPID,rownames(gtex_tpms_enz_log2))%>%.[!is.na(.)],] %>% as.matrix()

gtex_counts_no_zs <- apply(gtex_counts[,1:4], 2, function(x){
  x[which(is.infinite(x))] <- 0
  x})

gtex_counts2 <- as.data.frame(gtex_counts_no_zs) %>% data.frame(., tissue= new_gtex_annot$SMTS)

gtex_tissues <- levels(as.factor(gtex_counts2$tissue))

gtex_to_plot <- map(gtex_tissues, function(gts, g_df){
  df <- as.matrix(g_df[which(g_df$tissue==gts),-5])
  df <- apply(df,2,as.numeric)
  tpm_v_sel_med <- apply(df, 2, median)
  names(tpm_v_sel_med) <- colnames(df)
  tpm_v_sel_med
}, g_df=gtex_counts2) %>% do.call("rbind",.)

rownames(gtex_to_plot) <- gtex_tissues

pdf("../Results/Figures/bubbleheatmap_IDO1_TDO2_IDO2_IL4I1_GTEx.pdf",width = 8, height = 12, pointsize = 24)
corrplot(corr = t(as.matrix(gtex_to_plot)), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")
dev.off()