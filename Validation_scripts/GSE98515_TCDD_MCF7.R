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
tcdd100 <- read.delim("../Validation_datasets/GSE98515/GSE98515_raw_counts_100nM_TCDD.txt", sep = "\t")
tcdd10 <- read.delim("../Validation_datasets/GSE98515/GSE98515_counts_10nM_TCDD.txt", sep = "\t", stringsAsFactors = F)
raw_counts <- cbind(tcdd100, tcdd10[,grep("TCDD",colnames(tcdd10))]) %>% as.matrix()

## Metadata
metadatadf <- data.frame("sample" = c("C1", "C2", "C3", "C4","T1", "T2", "T3","T4", "d1", "d2", "d3", "d4") ,
                         "condition" = c("DMSO","DMSO", "DMSO", "DMSO", "TCDD", "TCDD", "TCDD", "TCDD","TCDD10", "TCDD10", "TCDD10", "TCDD10") )

rownames(metadatadf) <- make.names(metadatadf$sample, unique=TRUE)
metadatadf$sample <- NULL

## Create DGElist
dge_l <- DGEList(counts=raw_counts)
keep_l <- filterByExpr(dge_l)
dge_l <- dge_l[keep_l,,keep.lib.size=FALSE]
dge_l <- calcNormFactors(dge_l)
v_l <- voom(dge_l, plot=FALSE)

## Annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup2_2 <- getBM(mart=mart, attributes=c("ensembl_gene_id", "external_gene_name"),
                        filter="ensembl_gene_id", values=rownames(v_l$E), uniqueRows=TRUE)
rownames(v_l$E) <- annotLookup2_2$external_gene_name[match(rownames(v_l$E),annotLookup2_2$ensembl_gene_id)]

## Perform DGE
design <- model.matrix(~0+condition, data = metadatadf)
colnames(design) <- c("DMSO", "TCDD","TCDD10")

fit <- lmFit(v_l, design)
cm <- makeContrasts(TCDD-DMSO,TCDD10-DMSO, levels = design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2, robust = T)

tt_u_l <- topTable(fit2, coef=1, adjust="BH", number = Inf)

## Barcodeplot
index_vector_1 <- tt_u_l$ID %in% AHR_genes$Gene

pdf("../Results/Validations/GSE98515_MCF7/barcodeplot_MCF7.pdf",width = 12, height = 8)
barcodeplot(tt_u_l$t, index_vector_1, main = "GSE98514_TCDD_100nM_MCF7")
dev.off()

## ROAST
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
write.table(AHR_roast_res, "../Results/Validations/GSE98515_MCF7/GSE98515_roast_res.txt", sep = "\t")
