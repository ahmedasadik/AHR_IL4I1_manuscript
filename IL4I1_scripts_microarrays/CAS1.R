#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the analysis of the IL4I1 CRISPR/Cas9 KO in CAS1 cells.
#################################################
library(AffyGEx)

data("hgnc","msig.data.lists")

## add AHR signature genes
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", sep = "\t")

#### loading experimental condition covariates ####
exp_data_variables <- read.delim("../IL4I1_scripts_microarrays_metadata/CAS1_exp_data_variables.txt",sep = "\t")

setwd("../Results/IL4I1_microarrays/CAS1/")
dir.create("./RDS")
dir.create("./TopTables")
dir.create("./GSA")
dir.create("./Figures")

#### generate raw_cel dataset ####
raw_path <- "./CEL/"
list_cels <- list.files(raw_path, pattern = "CEL$", full.names = T)
raws_cel <- read.celfiles(list_cels)
saveRDS(raws_cel, "./RDS/raws_cel.rds")

#### Assigning sample names ####
s_names <- gsub("\\_\\HuGene\\-2\\_0\\-st\\_\\.CEL","",sampleNames(raws_cel)) %>%  gsub("GSM.([0-9]*_)","",.)
s_names <- gsub("CAS-1","CAS1", s_names) %>% gsub("-","_",.) %>%gsub("\\.","",.)
sampleNames(raws_cel) <- s_names

#### Create phenodata & metadata ####
condition <- gsub("CAS1_","",s_names)%>%gsub("(_.[0-9])","",.)
sample_id <- exp_data_variables$Sample
date_exp <- exp_data_variables$Date.performed
info <- data.frame(sample_id=sample_id, condition=condition, date_exp=date_exp)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'date_exp'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd

rma_2sd <- gen_rma_FUN(raws_cel, sd_val = 2,snames = s_names,ct_off = 0.75)

#### DGE ####
# Design matrix
design_eset <- model.matrix(~0+condition+date_exp, data=pd@data)
colnames(design_eset) <- gsub("condition","",colnames(design_eset))
contrast_matrix <- makeContrasts(IL4I1_KO_NTC=IL4I1_KO-NTC, IL4I1_KO_NTC01=IL4I1_KO-NTC01,NTC=NTC-NTC01,
                                 levels = design_eset)
## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_CAS1.rds")

## all tts
tts_2sd_u <- map(1:3, function(i){
  tt <- topTable(fit_2sd_eB, coef = i,  sort.by = "M",  number = Inf, adjust.method = "BH")
  tt <- tt[,c(29,25,40:45)]
})
names(tts_2sd_u) <- colnames(contrast_matrix)

## annotation update
library(dplyr)
tts_2sd_u[[1]] <- annot_hgcn_FUN(tts_2sd_u[[1]], hgnc = hgnc)
tts_2sd_u[[2]] <- annot_hgcn_FUN(tts_2sd_u[[2]], hgnc = hgnc)
tts_2sd_u[[3]] <- annot_hgcn_FUN(tts_2sd_u[[3]], hgnc = hgnc)

## add AHR signature genes
tts_2sd_u_AHR <- map(tts_2sd_u, function(x,y){
  x$AHR_target <- match(x$hGene, y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x
}, y=AHR_genes$Gene)

walk2(tts_2sd_u_AHR, paste("./TopTables/tt_2sd_",names(tts_2sd_u),"_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")

## Barcodeplots AHR signature
index_1 <- tts_2sd_u_AHR$IL4I1_KO_NTC$hGene%in%AHR_genes$Gene
index_2 <- tts_2sd_u_AHR$IL4I1_KO_NTC01$hGene%in%AHR_genes$Gene
index_3 <- tts_2sd_u_AHR$NTC$hGene%in%AHR_genes$Gene

pdf("./Figures/barcodeplot_CAS1.pdf")
barcodeplot(tts_2sd_u_AHR$IL4I1_KO_NTC$t,index_1, main="IL4I1_KO_vs_NCT")
barcodeplot(tts_2sd_u_AHR$IL4I1_KO_NTC01$t,index_2, main="IL4I1_KO_vs_NCT01")
barcodeplot(tts_2sd_u_AHR$NTC$t,index_3, main="NCT_vs_NCT01")
dev.off()

## Roast for the AHR signature
roast_FUN <- function(tt_ids, cont_mat,gs,rma_obj, des_mat, arr_wts=NULL){
  AHR_indexed_hi_lo <- ids2indices(gene.sets = gs, identifiers = tt_ids[,1])
  prbsets_gsa_2sd <- rownames(tt_ids)
  gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_obj@featureData@data$probesetid)
  gsa_2sd_eset <- rma_obj[gsa_2sd_idx]
  AHR_hi_lo_roast <- roast(oligo::exprs(gsa_2sd_eset),weights=arr_wts,index = AHR_indexed_hi_lo, design = des_mat,
                               contrast = cont_mat, set.statistic = "floormean")
  AHR_hi_lo_roast
}

all_roast <- map2(tts_2sd_u_AHR, as.data.frame(contrast_matrix), roast_FUN, gs=AHR_genes$Gene, rma_obj=rma_2sd, des_mat=design_eset) %>% do.call(rbind,.)
all_roast <- data.frame(Condition=rownames(all_roast), all_roast, stringsAsFactors = F)

write.table(all_roast, "./GSA/AHR_enrichment_CAS1_all.txt", row.names = F, sep = "\t")
