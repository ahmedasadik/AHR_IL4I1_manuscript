#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the analysis of the GSE102045 dataset.
#################################################

## Load libraries
library(purrr)
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(GEOquery)
library(org.Hs.eg.db)
library(biomaRt)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Load the raw data
GSE102045data <- read.delim("../Validation_datasets/GSE102045/GSE102045_raw_counts.csv", stringsAsFactors = F, sep = ",")

## Selecting samples
rownames(GSE102045data) <- GSE102045data$X
GSE102045data <- GSE102045data[,-1]
rem_ch22 <- grepl("CH22|Media", colnames(GSE102045data))
GSE102045data_1 <- GSE102045data[,!rem_ch22]

## Creating DGE list
dge_l <- DGEList(GSE102045data_1)

## Normalization
keep_l <- rowSums(dge_l$counts) >= 10
dge_l <- dge_l[keep_l,,keep.lib.size=FALSE]
dge_l <- calcNormFactors(dge_l)
v_l <- voom(dge_l, plot=FALSE)

## Annotations
mart <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup2_2 <- getBM(mart=mart, attributes=c("ensembl_gene_id", "external_gene_name"),
                        filter="ensembl_gene_id", values=rownames(v_l$E), uniqueRows=TRUE)
rownames(v_l$E) <- annotLookup2_2$external_gene_name[match(rownames(v_l$E),annotLookup2_2$ensembl_gene_id)]

## Performing DGE
donor <- factor(gsub("^JM..\\.","",colnames(v_l$E))%>%gsub("\\..*","",.))
group <- factor(gsub("^JM..\\.","",colnames(v_l$E))%>%
                  gsub("[A-Y]\\.","",.)%>%gsub("\\..*","",.))

design <- model.matrix(~0+group+donor)
colnames(design) <- gsub("group","",colnames(design))

fit <- lmFit(v_l, design)
cm <- makeContrasts(Th17_FICZ-Th17,levels = design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2, robust = T)
tt_u_l <- topTable(fit2, adjust="BH", number = Inf)

index_vector1 <- tt_u_l$ID %in% AHR_genes$Gene

## Barcodeplot
pdf("../Results/Validations/GSE102045_Th17/barcodeplot_Th17.pdf",width = 12, height = 8)
barcodeplot(tt_u_l$t, index_vector1, main="GSE102045_FICZ_TH17")
dev.off

## Roast
seq_roast_FUN <- function(v_eset, glist, des_mat, cont_mat, ...){
  index.vector <- rownames(v_eset$E) %in% glist
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

AHR_roast_res <- seq_roast_FUN(v_l, AHR_genes$Gene, design, cm[,1], robust=TRUE)
write.table(AHR_roast_res, "../Results/Validations/GSE102045_Th17/GSE102045_roast_res.txt", sep = "\t")
