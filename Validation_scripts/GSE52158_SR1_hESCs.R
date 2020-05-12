#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the analysis of the IL4I1 ectopic eexpression and IL4I1_shAHR U87MG cells.
#################################################

## Load libraries
library(purrr)
library(limma)
library(GEOquery)
library(Biobase)
library(R.utils)
library(hgu133plus2.db)
library(hgu133plus2cdf)
library(gcrma)
library(affy)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Metadata
GSE52158<- getGEO("GSE52158", GSEMatrix = FALSE)
gsm52158_metadata <- lapply(GSMList(GSE52158),function(x) {Meta(x)})
GSE52158_pt_dat <- as.data.frame(sapply(gsm52158_metadata, function(x)strsplit(x$title, ",")) %>% do.call("rbind",.))
timedesign <- GSE52158_pt_dat[-c(10:24),]
colnames(timedesign) <- c("group", "replicate")
timedesign$group <- factor(timedesign$group)
timedesign$replicate <- factor(timedesign$replicate)

## Read raw data
cel_raws <- ReadAffy(celfile.path = "../Validation_datasets/GSE52158/Raws/")
rma_cel <- expresso(cel_raws, bgcorrect.method = "rma",
                    normalize.method = "quantiles",
                    pmcorrect.method = "pmonly",
                    summary.method = "medianpolish")

sampleNames(rma_cel) <- rownames(timedesign)

edata <- exprs(rma_cel)

## Annotation
symbold <- mapIds(hgu133plus2.db, rownames(edata), "SYMBOL","PROBEID", multiVals = "list")
probes_mapped <- as.data.frame(unlist(map(symbold, paste, collapse=";")))
rownames(edata) <- probes_mapped[,1]

## Perform DGE
edata <- avereps(edata)

design <- model.matrix(~0+group, timedesign)
colnames(design) <- c("hESC", "SR1APS", "SR1DE")

fit <- lmFit(edata,design)
cont_mat <- makeContrasts(SR1APS-hESC,SR1DE-hESC, levels = design)
fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2, trend=TRUE, robust = TRUE)

tt1 <- topTable(fit2, coef = 1, number = Inf, sort.by = "M")

## Barcodeplot
index_vector_1 <- rownames(tt1) %in% AHR_genes$Gene
pdf("../Results/Validations/GSE52158_hESCs/barcodeplot_hESC.pdf", width = 12, height = 8)
barcodeplot(tt1$t, index_vector_1, main = "GSE52158_SR1APS_hESC")
dev.off()

## ROAST
affy_roast_FUN <- function(v_eset, glist, des_mat, cont_mat, ...){
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

GSE52158_SR1APS <- affy_roast_FUN(edata, AHR_genes$Gene, design, cont_mat[,1],robust=TRUE, trend=TRUE)
write.table(GSE52158_SR1APS, "../Results/Validations/GSE52158_hESCs/GSE52158_roast_res.txt", sep = "\t")
