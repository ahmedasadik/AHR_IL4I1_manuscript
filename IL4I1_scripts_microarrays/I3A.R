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
library(AffyGEx)

data("hgnc","msig.data.lists")

## add AHR signature genes
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", sep = "\t")

#### loading experimental condition covariates ####
exp_data_variables <- read.delim("../IL4I1_scripts_microarrays_metadata/I3A_exp_data_variables.txt", sep = "\t", stringsAsFactors = F)

setwd("../Results/IL4I1_microarrays/I3A/")
dir.create("./RDS")
dir.create("./TopTables")
dir.create("./GSA")
dir.create("./Figures")

## Raws
raw_path <- "./CEL/"
list_cels <- list.files(raw_path, pattern = "CEL$", full.names = T)
raws_cel <- read.celfiles(list_cels)
saveRDS(raws_cel,"./RDS/raws_cel.rds")

s_names <-gsub("\\_\\HuGene\\-2\\_0\\-st\\_\\.CEL","",sampleNames(raws_cel))  %>%  gsub("GSM.([0-9]*_)","",.)
sampleNames(raws_cel) <- s_names

#### Create phenodata & metadata ####
condition <- gsub("(\\U87_)","",s_names) %>% gsub("_0.","",.)
info <- data.frame(sample_id=s_names, condition=condition)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd
rma_2sd <- gen_rma_FUN(raws_cel, sd_val = 2, snames = s_names, ct_off = 0.25)

#### DGE ####
# Design matrix
design_eset <- model.matrix(~0+condition+exp_data_variables$Date_performed)
colnames(design_eset) <- gsub("condition","",colnames(design_eset))
colnames(design_eset) <- gsub("\\$","_", colnames(design_eset))
arr_wts <- arrayWeights(oligo::exprs(rma_2sd), design = design_eset)
contrast_matrix <- makeContrasts(I3A50uM-DMSO,levels = design_eset)

## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset, weights = arr_wts)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T, robust = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_U87_I3A.rds")

## all tts
tt_2sd_all <- topTable(fit_2sd_eB, number = Inf,adjust.method = "BH")
tt_2sd_u_all <- tt_2sd_all[,c(29,25,40:45)]

## Annotation update
require(dplyr)
tt_2sd_u_all <- annot_hgcn_FUN(tt_2sd_u_all, hgnc = hgnc)
saveRDS(tt_2sd_u_all, "./RDS/tt_2sd_u_annotated.rds")

## Assigning AHR genes
tt_2sd_u_AHR <- map(list(tt_2sd_u_all), function(x,y){
  x$AHR_target <- match(x[,1], y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x$logPV <- -log10(x$P.Value)
  x
}, y=AHR_genes$Gene)
walk2(list(tt_2sd_u_AHR), paste("./TopTables/tt_2sd_","I3A","_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")
names(tt_2sd_u_AHR) <- "I3A"

#### AHR signature enrichment
roast_FUN <- function(tt_ids, cont_mat, gs, rma_obj, des_mat, arr_wts,...){
  AHR_indexed_hi_lo <- ids2indices(gene.sets = gs, identifiers = tt_ids[,1])
  prbsets_gsa_2sd <- rownames(tt_ids)
  gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_obj@featureData@data$probesetid)
  gsa_2sd_eset <- rma_obj[gsa_2sd_idx]
  AHR_hi_lo_roast <- roast(oligo::exprs(gsa_2sd_eset),weights=arr_wts,index = AHR_indexed_hi_lo, design = des_mat,
                               contrast = cont_mat, set.statistic = "floormean")
  AHR_hi_lo_roast
}

all_roast <- map2(tt_2sd_u_AHR, as.data.frame(contrast_matrix), roast_FUN, gs=AHR_genes$Gene, rma_obj=rma_2sd, des_mat=design_eset,arr_wts=arr_wts, trend=T, robust=T) %>% do.call(rbind,.)
all_roast <- data.frame(Condition=rownames(all_roast), all_roast, stringsAsFactors = F)
write.table(all_roast, "./GSA/AHR_enrichment_I3A.txt", row.names = F, sep = "\t")

pdf("./Figures/barcodeplot_I3A_U87.pdf", width = 12, height = 8)
index_vector <- tt_2sd_u_AHR$I3A$hGene %in% AHR_genes$Gene
barcodeplot(tt_2sd_u_AHR$I3A$t, index_vector, main=all_roast$Condition[1])
dev.off()
