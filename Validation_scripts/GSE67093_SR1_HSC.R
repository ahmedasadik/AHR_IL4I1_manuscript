#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the analysis of the IL4I1 ectopic eexpression and IL4I1_shAHR U87MG cells.
#################################################

## Load libraries
library(GEOquery)
library(limma)
library(dplyr)
library(purrr)
library(factoextra)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Load raw data
all_raws2 <- read.delim("../Validation_datasets/GSE67093/GSE67093_non-normalized.txt",
                        sep = "\t", skip= 8, stringsAsFactors = F, header = TRUE)

## Metadata
GSE67093 <- getGEO("GSE67093", GSEMatrix = FALSE)
gsm_metadata <- lapply(GSMList(GSE67093),function(x) {Meta(x)})
GSE67093_pt_dat <- map(gsm_metadata, function(x){
  metdat <- strsplit(x$title, split = "-")[[1]]
  if (length(metdat) == 2){
    metdat <- metdat[c(1,2)]
  } else {
    metdat <- metdat[c(1,4)]
  }
  metdat
}
) %>% do.call("rbind",.)

GSE67093_pt_dat <- as.data.frame(GSE67093_pt_dat, stringsAsFactors = FALSE)
GSE67093_pt_dat <- GSE67093_pt_dat[!GSE67093_pt_dat$V1 == "Pre",]
use_rownames <- rownames(GSE67093_pt_dat) 
GSE67093_pt_dat <- strsplit(GSE67093_pt_dat$V1, " ")  %>% do.call("rbind",.)
GSE67093_pt_dat <- as.data.frame(GSE67093_pt_dat, stringsAsFactors = FALSE)
rownames(GSE67093_pt_dat) <- use_rownames
GSE67093_pt_dat$V1 <- paste0(GSE67093_pt_dat$V1, GSE67093_pt_dat$V2, GSE67093_pt_dat$V3, GSE67093_pt_dat$V4) 
GSE67093_pt_dat <- GSE67093_pt_dat[,-c(2:4)]
GSE67093_pt_dat$V1 <- gsub(",", "_", GSE67093_pt_dat$V1)
GSE67093_pt_dat[,2] <- paste("HT",c(2:5,7:10,14:17,19:22), sep = "_")
GSE67093_pt_dat[,3] <- rep(paste("donor",1:4,sep = ""),each=4)
colnames(GSE67093_pt_dat) <- c("condition", "sample", "donor")

## Create EListRaw
all_raws2 <- all_raws2[, -grep("^HT.1.1", colnames(all_raws2))]
all_raws2 <- all_raws2[, -grep("^HT18", colnames(all_raws2))]
all_raws2 <- all_raws2[, -grep("^HT13", colnames(all_raws2))]
all_raws2 <- all_raws2[, -grep("^HT_06", colnames(all_raws2))]
all_raws2 <- all_raws2[, -grep("^HT_\\d1", colnames(all_raws2))]
all_raws2 <- all_raws2[, -grep("^HT_12", colnames(all_raws2))]

obj_arrays <- new("EListRaw")
obj_arrays$source <- "illumina"
obj_arrays$E <- all_raws2[,grep("AVG_Signal",colnames(all_raws2))] %>% as.matrix()
obj_arrays$gene <- data.frame(SYMBOL=all_raws2$SYMBOL, stringsAsFactors = F)
obj_arrays$other <- list()
obj_arrays$other$Detection <- all_raws2[,grep("Detection",colnames(all_raws2))] %>% as.matrix()
rownames(obj_arrays$E) <- rownames(obj_arrays$other$Detection) <- obj_arrays$gene$SYMBOL
arr_norm <- neqc(obj_arrays) 
obj_arrays2 <- arr_norm
good_probes <- rowSums(obj_arrays2$other$Detection) <= 0.05
obj_arrays2$E <- obj_arrays2$E[good_probes,]
obj_arrays2$other$Detection <-  obj_arrays2$other$Detection[good_probes,]
obj_arrays2$gene <- obj_arrays2$gene[good_probes,]

## Perform DGE
design <- model.matrix(~0+condition+donor, GSE67093_pt_dat)
colnames(design) <- gsub("condition","",colnames(design)) %>% gsub("^donor","",.)%>%gsub(" ","",.)

fit <- lmFit(obj_arrays2,design)
cont_mat <- makeContrasts("withSR1_withoutZFN - withoutSR1_withoutZFN",
                          "withoutSR1_withZFN - withoutSR1_withoutZFN",
                          "withSR1_withZFN - withoutSR1_withoutZFN",levels = design)
fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2, trend=TRUE, robust = TRUE)
tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "M")

## Barcodeplot
index_vector_1 <- tt$ID %in% AHR_genes$Gene
pdf("../Results/Validations/GSE67093_HSC/barcodeplot_HSC.pdf",width = 12, height = 8)
barcodeplot(tt$t, index_vector_1, main="GSE67093_SR1_HSC")
dev.off()

## ROAST
ilmn_roast_FUN <- function(v_eset, glist, des_mat, cont_mat, ...){
  index.vector <- rownames(v_eset) %in% glist
  AHR_roast <- roast(v_eset, index = index.vector, design = des_mat,
                     contrast = cont_mat, set.statistic = "floormean")
  roast_cols <- c("NGenes",	"PropDown",	"PropUp",	"Direction",	"PValue",	"FDR",	"PValue.Mixed",	"FDR.Mixed")
  AHR_roast_df <- as.matrix(rep(0,length(roast_cols)),nrow=1) %>% t() %>% as.data.frame()
  colnames(AHR_roast_df) <- roast_cols
  AHR_roast_df$NGenes <- AHR_roast$ngenes.in.set
  AHR_roast_df$PropDown <- AHR_roast$p.value$Active.Prop[1] %>% round(.,4)
  AHR_roast_df$PropUp <- AHR_roast$p.value$Active.Prop[2] %>% round(.,4)
  AHR_roast_df$Direction <- ifelse(AHR_roast_df$PropDown > AHR_roast_df$PropUp, "Down", "Up")
  AHR_roast_df$PValue <- AHR_roast_df$FDR <- ifelse(AHR_roast_df$Direction=="Down",
                                                    round(AHR_roast$p.value$P.Value[3], 4),
                                                    round(AHR_roast$p.value$P.Value[3],4))
  AHR_roast_df$PValue.Mixed <- AHR_roast_df$FDR.Mixed <- round(AHR_roast$p.value$P.Value[4],4)
  AHR_roast_df 
}

GSE67093_ROAST_SR1_vs_no_SR1_all_noZNF <- ilmn_roast_FUN(obj_arrays2$E, AHR_genes$Gene, design, cont_mat = cont_mat[,1],set.statistic="floormean")
write.table(GSE67093_ROAST_SR1_vs_no_SR1_all_noZNF, "../Results/Validations/GSE67093_HSC/GSE67093_roast_res.txt", sep = "\t")
