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
library(dplyr)
library(limma)
library(purrr)
library(biomaRt)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Metadata
GSE32026 <- getGEO("GSE32026", GSEMatrix = FALSE)
gsm32026_metadata <- lapply(GSMList(GSE32026),function(x) {Meta(x)})
pt_dat_32026 <- as.data.frame(sapply(gsm32026_metadata, function(x)strsplit(x$title, " ")) %>% do.call("rbind",.))
remove <- c(2, 3, 5, 6)
pt_dat_32026 <- pt_dat_32026[,-remove]
pt_dat_32026$V7 <- paste0(pt_dat_32026$V7, pt_dat_32026$V8)
pt_dat_32026 <- pt_dat_32026[,-4]
colnames(pt_dat_32026) <- c("condition", "type", "replicate")
pt_dat_32026$condition <- gsub("-", "", pt_dat_32026$condition)
pt_dat_32026 <- pt_dat_32026[-c(5:12),]
pt_dat_32026 <- pt_dat_32026[-c(9:16),]
pt_dat_32026[5:8,1] <- paste0(pt_dat_32026[5:8,1], "DIFF")

## Load targets file
targets <- read.delim("../Validation_datasets/GSE32026/targets.txt", sep = "\t", stringsAsFactors = F)

## Read raw files
RG <- read.maimages(targets, path = "../Validation_datasets/GSE32026/Raws/", source = "genepix")
y <- backgroundCorrect(RG, method="normexp", offset=50)
y <- normalizeWithinArrays(y, method="loess")

## Annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
results<- getBM(attributes = c("refseq_mrna","hgnc_symbol"), filters="refseq_mrna", values=RG$genes$Name, mart=mart)
rownames(y) <- y$genes$Name
y$genes$symbol <- results$hgnc_symbol[match(y$genes$Name, results$refseq_mrna)]
y$genes <- na.omit(y$genes)
y$M <- y$M[y$genes$Name,]
y$A <- y$A[y$genes$Name,]
rownames(y) <- y$genes$symbol

## Filtering
IsExprM <- rowSums(y$M > 0) >= 4
yfilt <- y[IsExprM,]
y.ave <- avereps(yfilt)

## DGE
design <- model.matrix(~0+condition, pt_dat_32026)
colnames(design) <- c("TCDD", "TCDD_DIFF")
fit <- lmFit(y.ave, design)
fit2 <- eBayes(fit)
tt1 <- topTable(fit2, coef = 1,  adjust="BH", number = Inf)

## Barcodeplots
index_vector_1 <- tt1$symbol %in% AHR_genes$Gene
pdf("../Results/Validations/GSE32026_hMADs/barcodeplot_hMADs.pdf",width = 12, height = 8)
barcodeplot(tt1$t, index_vector_1, main="GSE32026_TCDD_hMADS")
dev.off()

## ROAST
agilent_roast_FUN <- function(v_eset, glist, des_mat, cont_mat, ...){
  index.vector <- v_eset$genes$symbol %in% glist
  AHR_roast <- roast(v_eset, index = index.vector, design = des_mat,
                     contrast = cont_mat, set.statistic = "floormean")
  roast_cols <- c("NGenes",	"PropDown",	"PropUp",	"Direction",	"PValue",	"FDR",	"PValue.Mixed",	"FDR.Mixed")
  AHR_roast_df <- as.matrix(rep(0,length(roast_cols)),nrow=1) %>% t() %>% as.data.frame()
  colnames(AHR_roast_df) <- roast_cols
  AHR_roast_df$NGenes <- AHR_roast$ngenes.in.set
  AHR_roast_df$PropDown <- AHR_roast$p.value$Active.Prop[1] %>% round(.,5)
  AHR_roast_df$PropUp <- AHR_roast$p.value$Active.Prop[2] %>% round(.,5)
  AHR_roast_df$Direction <- ifelse(AHR_roast_df$PropDown > AHR_roast_df$PropUp, "Down", "Up")
  AHR_roast_df$PValue <- AHR_roast_df$FDR <- ifelse(AHR_roast_df$Direction=="Down",
                                                    round(AHR_roast$p.value$P.Value[3], 5),
                                                    round(AHR_roast$p.value$P.Value[3],5))
  AHR_roast_df$PValue.Mixed <- AHR_roast_df$FDR.Mixed <- round(AHR_roast$p.value$P.Value[4],5)
  AHR_roast_df 
}

GSE32026_ROAST <- agilent_roast_FUN(y.ave, AHR_genes$Gene, design, cont_mat = 1)
write.table(GSE32026_ROAST,"../Results/Validations/GSE32026_hMADs/GSE32026_roast_res.txt", sep = "\t")
