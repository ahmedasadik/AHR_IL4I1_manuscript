#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the analysis of the relapsed CLL patient dataset GSE58211
#################################################

## Load libraries
library(GEOquery)
library(purrr)
library(AffyGEx)
library(huex10stprobeset.db)
library(R.utils)
library(ComplexHeatmap)
library(circlize)

## Source the functions and parameters files
source("../functions_and_parameters.R")

## Create directory
dir.create("../Results/CLL_GSE58211")

# dir.create("../Results/CLL_GSE58211/CEL")
# dir.create("../Results/CLL_GSE58211/RDS")
dir.create("../Results/CLL_GSE58211/TopTables")
# dir.create("../Results/CLL_GSE58211/GSA")
dir.create("../Results/CLL_GSE58211/Figures")
# dir.create("../Results/CLL_GSE58211/TFSDB")

# ## Downloading the raw files
# GSE58211_files <- getGEOSuppFiles("GSE58211",baseDir = "../Results/CLL_GSE58211")
# untar("../Results/CLL_GSE58211/GSE58211/GSE58211_RAW.tar", tarfile = ,exdir = "../Results/CLL_GSE58211/CEL")
# tar_list <- list.files("../Results/CLL_GSE58211/CEL", pattern = "GSM", full.names = T)
# map(tar_list,gunzip)
# 
# ## Metadata
# GSE58211 <- getGEO("GSE58211", GSEMatrix = FALSE)
# GSE58211_metadata <- lapply(GSMList(GSE58211),function(x) {Meta(x)})
# GSE58211_pt_dat <- as.data.frame(sapply(GSE58211_metadata, function(x)strsplit(x$title, " ")) %>% do.call("rbind",.))
# GSE58211_pt_dat$V3 <- paste(GSE58211_pt_dat$V3,GSE58211_pt_dat$V6,sep = "_")%>%gsub(" |\\#","",.)
# GSE58211_pt_dat <- GSE58211_pt_dat[,-c(4:6)]
# colnames(GSE58211_pt_dat) <- c("condition","disease","pt_ID")
# GSE58211_pt_dat$GSM <- rownames(GSE58211_pt_dat)
# GSE58211_pt_dat$binet <- sapply(GSE58211_metadata, function(x){x$characteristics_ch1[3]%>%gsub("\\:","",.)%>%gsub(" ","_",.)})
# 
# cel_files <- list.files("../Results/CLL_GSE58211/CEL/", pattern = "\\.CEL$", full.names = T)
# raw_cel <- read.celfiles(cel_files)
# 
# saveRDS(raw_cel,"../Results/CLL_GSE58211/RDS/CLL_GSE58211_raw_cel.rds")
# 
# ## Assigning sample names
# s_names <- gsub(".//","",sampleNames(raw_cel))%>%gsub("_.*","",.)
# sampleNames(raw_cel) <- s_names
# 
# #### Create phenodata & metadata ####
# condition <- GSE58211_pt_dat$condition
# sample_id <- GSE58211_pt_dat$GSM
# donor <- GSE58211_pt_dat$pt_ID
# info <- data.frame(sample_id=sample_id, condition=condition, donor=donor)
# rownames(info) <- s_names
# metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'donor'), channel=factor('_ALL_'))
# pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
# phenoData(raw_cel) <- pd
# 
# rma_2sd <- gen_rma_FUN(raw_cel, sd_val = 0.5, snames = s_names, ct_off = 0.999, rds_path ="../Results/CLL_GSE58211/RDS/CLL_GSE58211_rma_2sd.rds")

## load normalized expression values
rma_2sd <- readRDS("../Resources/GSE58211/CLL_GSE58211_rma_2sd.rds")

## Estimate the AHR score for the CLL patients
AHR_signature <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

rma_exprs <- exprs(rma_2sd)
rownames(rma_exprs) <- rma_2sd@featureData@data$geneid

# ## Generating an AHR score
# rma_gsva <- GSVA::gsva(rma_exprs,list(AHR_signature$Gene),parallel.sz=no_cores)
# saveRDS(rma_gsva,"../Results/CLL_GSE58211/RDS/gsva_AHR.rds")

rma_gsva <- readRDS("../Resources/GSE58211/CLL_gsva_AHR.rds")

df <- data.frame(sample_id=colnames(rma_gsva),AHR_score=rma_gsva[1,],IL4I1=rma_exprs[which(rownames(rma_exprs)=="IL4I1"),], stringsAsFactors = F)
df2 <- dplyr::full_join(df,rma_2sd@phenoData@data,"sample_id")
df2_melt <- reshape::melt(df2)

# ## Run WGCNA
# library(WGCNA)
# 
# # Enabling WGCNA threading
# enableWGCNAThreads(nThreads = no_cores)
# 
# # Estimating soft threshold power
# # Selecting the powers for soft thresholding
# powers <- c(c(1:10), seq(from = 12, to=40, by=2))
# 
# # Select soft threshold
# CLL_sft <- pickSoftThreshold(t(rma_exprs), dataIsExpr = TRUE, networkType = "signed hybrid",
#                     blockSize = dim(t(rma_exprs))[[2]], powerVector = powers, verbose = 5)
# saveRDS(CLL_sft, "../Results/CLL_GSE58211/RDS/CLL_sft.rds")
# 
# ## RUN WGCNA and generate coexpression modules
# CLL_TOMs <- blockwiseModules(datExpr = t(rma_exprs), corType = "bicor",
#                    maxBlockSize = dim(t(rma_exprs))[[2]], nThreads = 30,
#                    randomSeed = 0861, power = CLL_sft$powerEstimate,
#                    networkType = "signed hybrid",
#                    TOMType = "signed", TOMDenom = "mean", saveTOMs = TRUE,
#                    saveTOMFileBase = "GSE58211",
#                    pamStage = TRUE, pamRespectsDendro = FALSE,verbose = 3)
# saveRDS(CLL_TOMS_new,"../Results/CLL_GSE58211/RDS/CLL_GSE58211_TOMs.rds")
# 
# ## Extract the gene names for the different modules
# CLL_modules_genes <- data.frame(genes=rownames(rma_exprs), module=CLL_TOMs_new$colors, stringsAsFactors = F)
# saveRDS(CLL_modules_genes, "../Results/CLL_GSE58211/RDS/CLL_module_genes.rds")

CLL_TOMs <- readRDS("../Resources/GSE58211/CLL_GSE58211_TOMs.rds")
CLL_modules_genes <- readRDS("../Resources/GSE58211/CLL_module_genes.rds")

## Correlating between the modules and AHR activation?
cor_mat <- Hmisc::rcorr(as.matrix(CLL_TOMs$MEs), as.numeric(rma_gsva[1,]))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    modules = rownames(cormat)[row(cormat)[ut]],
    cor_genes = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut], stringsAsFactors = FALSE)
}
cor_df <- flattenCorrMatrix(cor_mat$r, cor_mat$P)
cor_df$cor_genes[which(cor_df$cor_genes == "y")] <- "AHR_score"
cor_df <- cor_df[cor_df$cor_genes=="AHR_score",] %>%
  .[which(.$cor>=0.2 | .$cor<=(0.2*-1)),] %>%
  .[.$p <=0.05,] %>% .[order(.$cor, decreasing = TRUE),]
cor_df$modules <- gsub("^ME","",cor_df$modules)
cor_df <- data.frame(cor_df, stringsAsFactors = F)
cor_df$Has_IL4I1 <- "no"
cor_df$Has_IL4I1[which(cor_df$modules=="grey")]<-"yes"
cor_df <- cor_df[order(cor_df$cor,decreasing = T),]
write.csv(cor_df,"../Results/CLL_GSE58211/TopTables/CLL_WGCNA_modules_cor_AHR.csv")

pdf("../Results/CLL_GSE58211/Figures/CLL_AHR_IL4I1_correlation.pdf",width = 12, height = 8)
ggscatter(df2,"AHR_score","IL4I1",add="reg.line",cor.coef = T)
ggbarplot(cor_df,y = "cor",x="modules",fill = "Has_IL4I1",palette = c("black","red"))
dev.off()

## Generating a circos plot
ordered_cols <- c("grey","blue","pink","red")
factors2 <- as.factor(rep(c("IL4I1", "IDO1","IDO2","TDO2","DDC","TPH1","TPH2", ordered_cols), 
                          c(100,100,100,100,100,100,100, unlist(map(ordered_cols,function(x){
                            length(CLL_modules_genes$module[CLL_modules_genes$module==x])})))
))

factors2 <- factor(factors2, levels = c("grey","blue","red","pink","IL4I1","IDO1","IDO2","TDO2","DDC","TPH1","TPH2"))
to_subset <- as.character(levels(factors2))

pdf("../Results/Figures/CLL_WGCNA_circos_plot_7.pdf", width = 8, height = 8)
## plot outline
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 4)
x_limits <- cbind(rep(0,length(to_subset)), table(factors2))
circos.initialize(to_subset, xlim = x_limits)
circos.track(ylim = c(0, 0.5), bg.border = NA, track.height=0.2,
             panel.fun = function(x, y) {
               x_p <- get.cell.meta.data("xlim")
               y_p <- get.cell.meta.data("ylim")
               circos.rect(x_p[1], y_p[1],x_p[2], y_p[2],
                           col = ifelse(CELL_META$sector.index %in% colors(),CELL_META$sector.index, "white"),
                           border = ifelse(CELL_META$sector.index %in% colors(),"black", NA))
             })
## add enzyme names
circos.text(50, 0, labels="IL4I1", sector.index = "IL4I1", facing = "clockwise", niceFacing = T, col = "black", adj = c(-0.05, degree(1)))
circos.text(50, 0, labels="IDO1", sector.index = "IDO1", facing = "clockwise", niceFacing = T, col = "black", adj = c(-0.05, degree(1)))
circos.text(50, 0, labels="IDO2", sector.index = "IDO2", facing = "clockwise", niceFacing = T, col = "black", adj = c(-0.05, degree(1)))
circos.text(50, 0, labels="TDO2", sector.index = "TDO2", facing = "clockwise", niceFacing = T, col = "black", adj = c(-0.05, degree(1)))
circos.text(50, 0, labels="DDC", sector.index = "DDC", facing = "clockwise", niceFacing = T, col = "black", adj = c(-0.05, degree(1)))
circos.text(50, 0, labels="TPH1", sector.index = "TPH1", facing = "clockwise", niceFacing = T, col = "black", adj = c(-0.05, degree(1)))
circos.text(50, 0, labels="TPH2", sector.index = "TPH2", facing = "clockwise", niceFacing = T, col = "black", adj = c(-0.05, degree(1)))
circos.link("IL4I1",point1 = c(0,100), "grey", point2 =c(12880,12980),  col = "grey",border = "black")
circos.link("IDO1",point1 = c(0,100), "grey", point2 =c(12880,12980),  col = "grey",border = "black")
circos.link("IDO2",point1 = c(0,100), "grey", point2 =c(12880,12980),  col = "grey",border = "black")
circos.link("TDO2",point1 = c(0,100), "grey", point2 =c(12880,12980),  col = "grey",border = "black")
circos.link("DDC",point1 = c(0,100), "grey", point2 =c(12880,12980),  col = "grey",border = "black")
circos.link("TPH1",point1 = c(0,100), "grey", point2 =c(12880,12980),  col = "grey",border = "black")
circos.link("TPH2",point1 = c(0,100), "grey", point2 =c(12880,12980),  col = "grey",border = "black")
title("CLL_GSE58211",cex=0.5)
circos.clear()
dev.off()
