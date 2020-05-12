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
library(gcrma)
library(affy)
library(stringr)
library(dplyr)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Metadata
GSE28878 <- getGEO("GSE28878", GSEMatrix = FALSE)
gsm28878_metadata <- lapply(GSMList(GSE28878),function(x) {Meta(x)})
GSE28878_pt_dat <- as.data.frame(sapply(gsm28878_metadata, function(x)strsplit(x$title, " ")) %>% do.call("rbind",.))
remove <- c(1, 3:5, 7, 9:10)
GSE28878_pt_dat <- GSE28878_pt_dat[,-remove]
GSE28878_pt_dat <-  subset(GSE28878_pt_dat, V6 == "BaP" | V6 == "TCDD" | V6 == "DMSO" & V2 == "A,")
GSE28878_pt_dat$V8 <- gsub("," , "", GSE28878_pt_dat$V8)
GSE28878_pt_dat$V2 <- gsub("," , "", GSE28878_pt_dat$V2)
colnames(GSE28878_pt_dat) <- c("Serie", "condition", "timepoint", "replicate")
GSE28878_pt_dat$condition <-gsub("DMSO", "ADMSO"  , GSE28878_pt_dat$condition)
GSE28878_pt_dat$batch <- rep(paste("b", 1:3, sep = "") %>% as.factor())

## Read raw files
cel_raws <- ReadAffy(celfile.path = "../Validation_datasets/GSE28878/Raws/")
rma_cel <- expresso(cel_raws, bgcorrect.method = "rma",
                    normalize.method = "quantiles",
                    pmcorrect.method = "pmonly",
                    summary.method = "medianpolish")

sampleNames(rma_cel) <- rownames(GSE28878_pt_dat)
edata <- exprs(rma_cel)

## Annotation
symbold <- mapIds(hgu133plus2.db, rownames(edata), "SYMBOL","PROBEID", multiVals = "list")
probes_mapped <- data.frame(unlist(map(symbold, paste, collapse=";")), stringsAsFactors = F)
rownames(edata) <- probes_mapped[,1]

## Averaging expression matrix
edata_av <- avereps(edata)

## DGE
design <- model.matrix(~0+condition:timepoint, GSE28878_pt_dat)
colnames(design) <- c("DMSO_12", "BaP_12", "TCDD_12",
                      "DMSO_24", "BaP_24", "TCDD_24",
                      "DMSO_48", "BaP_48", "TCDD_48")

fit <- lmFit(edata_av, design)
cont.wt <- makeContrasts("TCDD_12-DMSO_12", "BaP_12-DMSO_12", "TCDD_24-DMSO_24",
                         "BaP_24-DMSO_24", "TCDD_48-DMSO_48", "BaP_48-DMSO_48", levels=design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)

tt3 <-  topTable(fit2, coef=3, adjust="BH", number = Inf)

## Barcodeplot
index_vector_3 <- rownames(tt3) %in% AHR_genes$Gene
pdf("../Results/Validations/GSE28878_HepG2/barcodeplot_HepG2.pdf",width = 12, height = 8)
barcodeplot(tt3$t, index_vector_3, main= "GSE28878_TCDD_HepG2_24H")
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

GSE28878_ROAST_TCDD24 <- affy_roast_FUN(edata_av, AHR_genes$Gene, design, cont.wt[,3])
write.table(GSE28878_ROAST_TCDD24,"../Results/Validations/GSE28878_HepG2/GSE28878_roast_res.txt", sep = "\t")
