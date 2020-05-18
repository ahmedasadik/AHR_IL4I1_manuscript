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
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(GEOquery)
library(org.Hs.eg.db)
library(biomaRt)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Load raw data
raw_paths <- list.files("../Validation_datasets/GSE48843/Raws/", full.names = T)
raw_names <- list.files("../Validation_datasets/GSE48843/Raws/")
raws_read <- map(raw_paths,function(x){
  df <- read.delim(x,stringsAsFactors=F,header=F)
  colnames(df) <- c("Chr","Start", "Stop","Gene","RPKM","Count")
  data.frame(df[,c("Gene","Count")],stringsAsFactors = F)
  })
names(raws_read) <- gsub("_.*","",raw_names)

counts_df <- map(raws_read,function(x)x[,2])%>%do.call(cbind,.)
rownames(counts_df) <- raws_read$GSM1202511$Gene

## Metadata
GSE48843<- getGEO("GSE48843", GSEMatrix = FALSE)
GSE48843_metadata <- lapply(GSMList(GSE48843),function(x) {Meta(x)})
GSE48843_pt_dat <- as.data.frame(sapply(GSE48843_metadata, function(x){
  v1 <- strsplit(x$characteristics_ch1, ": ") %>% map(.,function(a){a[[2]]}) %>% unlist()
  c(v1,x$title)
}))%>%t() %>% as.data.frame()
colnames(GSE48843_pt_dat) <- c("treatments","time","source","sample_ID")
GSE48843_pt_dat$donor <- gsub("^([0-9]-)|^([0-9].-)","",GSE48843_pt_dat$sample_ID)%>%gsub("\\-.*","",.)
GSE48843_pt_dat2 <- GSE48843_pt_dat[grep("^DMSO|^(SR1_500nM)",GSE48843_pt_dat$treatments),] %>% .[-grep("BIO",.$treatments),]
GSE48843_pt_dat2$group <- ifelse(grepl("DMSO",GSE48843_pt_dat2$treatments),"DMSO","SR1")

## Create DGE
counts_df2 <- counts_df[,match(rownames(GSE48843_pt_dat2),colnames(counts_df))]
dge_l <- DGEList((counts_df2))
keep_l <- filterByExpr(dge_l)
dge_l <- dge_l[keep_l,,keep.lib.size=FALSE]
dge_l <- calcNormFactors(dge_l)
v_l <- voom(dge_l, plot=FALSE)

## Perform DGE
design <- model.matrix(~0+group+donor, GSE48843_pt_dat2)
colnames(design) <- gsub("group","",colnames(design)) %>% gsub("\\(|\\)","",.)

fit <- lmFit(v_l, design)
cm <- makeContrasts(SR1-DMSO, levels = design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)
tt_u_l <- topTable(fit2, coef = 1, adjust="BH", number = Inf)

## Barcodeplot
index_vector1 <- tt_u_l$ID%in% AHR_genes$Gene
pdf("../Results/Validations/GSE48843_AML/barcodeplot_AML.pdf",width = 12, height = 8)
barcodeplot(tt_u_l$t, index_vector1, main="GSE48843_SR1_AML")
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
AHR_roast_res <- seq_roast_FUN(v_l, AHR_genes$Gene, design, cm[,1])

write.table(AHR_roast_res, "../Results/Validations/GSE48843_AML/GSE102045_roast_res.txt", sep = "\t")
