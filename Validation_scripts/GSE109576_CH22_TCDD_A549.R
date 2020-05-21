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
library(illuminaHumanv4.db)
library(GEOquery)
library(dplyr)
library(limma)
library(purrr)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Load the raw data
allraws <- read.delim("../Validation_datasets/GSE109576/GSE109576_non_normalized.txt", sep = "\t", stringsAsFactors = F, header = TRUE)

#Metadata
GSE109576 <- getGEO("GSE109576", GSEMatrix = FALSE)

#Patient Metadata
gsm_109576_metadata <- lapply(GSMList(GSE109576),function(x) {Meta(x)})

GSE109576_pt_dat <- as.data.frame(sapply(gsm_109576_metadata, function(x)strsplit(x$title, " ")) %>% do.call("rbind",.))

remove <- c(1:4, 6 )

GSE109576_pt_dat<- GSE109576_pt_dat[,-remove]

GSE109576_pt_dat$V7 <- gsub(".*_", " ", GSE109576_pt_dat$V7)

colnames(GSE109576_pt_dat) <- c("Condition", "Replicate")

GSE109576_pt_dat$Condition <-   gsub("DMSO", "Control ", GSE109576_pt_dat$Condition)

## Preparing the raw data
allraws <- allraws[,-c(26:30)]
rownames(allraws) <- allraws$ID_REF
allraws <- allraws[,-1]
symbold <- mapIds(illuminaHumanv4.db, rownames(allraws), "SYMBOL","PROBEID", multiVals = "list")
probes_mapped <- data.frame(unlist(map(symbold, paste, collapse=";")), stringsAsFactors = F)

## Create EListRaw
obj_arrays <- new("EListRaw")
obj_arrays$source <- "illumina"
obj_arrays$E <- allraws[,seq(1, ncol(allraws), by = 2)] %>% as.matrix()
colnames(obj_arrays$E) <- rownames(GSE109576_pt_dat)
obj_arrays$gene <- data.frame(SYMBOL=probes_mapped[,1], stringsAsFactors = F)
rownames(obj_arrays$gene) <- rownames(probes_mapped)
obj_arrays$other <- list()
obj_arrays$other$Detection <- allraws[,seq(2, ncol(allraws), by = 2)] %>% as.matrix()
colnames(obj_arrays$other$Detection) <- colnames(obj_arrays$E)

## Normalization
arr_norm <- neqc(obj_arrays, offset = 50)

## remove un-annotated probes 
na_index <- grepl("NA", arr_norm$gene$SYMBOL)

arr_norm2 <- arr_norm[!na_index,]
arr_norm2$gene <- arr_norm2$gene[!na_index,]

good_probes <- rowSums(arr_norm2$other$Detection < 0.05) >= 3

arr_norm2$E <- arr_norm2$E[good_probes,]
arr_norm2$gene <- arr_norm2$gene[good_probes]
rownames(arr_norm2$E) <- arr_norm2$gene
arr_norm2$other$Detection <- arr_norm2$other$Detection[good_probes,]

design <- model.matrix(~0+Condition, GSE109576_pt_dat)

colnames(design) <- gsub("Condition","",colnames(design)) %>% gsub("\\-","",.)%>%gsub(" ","",.)

fit <- lmFit(arr_norm2, design)

cont_mat <- makeContrasts(CH223191-Control, TCDD-Control, TCDDCH223191-Control, levels = design)

fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2, trend=TRUE, robust = TRUE)

tt1 <- topTable(fit2, coef = 1, number = Inf, sort.by = "M")
tt2 <- topTable(fit2, coef = 2, number = Inf, sort.by = "M")
tt3 <- topTable(fit2, coef = 3, number = Inf, sort.by = "M")


## Barcode plots
index_vector_1 <- tt1$ID %in% AHR_genes$Gene
index_vector_2 <- tt2$ID %in% AHR_genes$Gene

pdf("../Results/Validations/GSE109576_A549/barcodeplots_A549.pdf",width = 12, height = 8)
barcodeplot(tt1$t, index_vector_1, main = "GSE109576_CH223139_A549")
barcodeplot(tt2$t, index_vector_2, main = "GSE109576_TCDD_A549")
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

cont1_roast <- ilmn_roast_FUN(arr_norm2$E, AHR_genes$Gene, design, cont_mat[,1], trend=TRUE, robust=TRUE)
cont2_roast <- ilmn_roast_FUN(arr_norm2$E, AHR_genes$Gene, design, cont_mat[,2], trend=TRUE, robust=TRUE)

df <- rbind(cont1_roast, cont2_roast)
df <- data.frame(Condition=c("CH22","TCDD"),df,stringsAsFactors = F)

write.table(df,"../Results/Validations/GSE109576_A549/GSE109576_roast_res.txt")
  