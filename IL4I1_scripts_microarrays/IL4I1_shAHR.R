#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the analysis of the IL4I1 ectopic expression and IL4I1_shAHR U87MG cells.
#################################################
library(AffyGEx)

data("hgnc","msig.data.lists")

## add AHR signature genes
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", sep = "\t")

#### loading experimental condition covariates ####
exp_data_variables <- read.delim("../IL4I1_scripts_microarrays_metadata/IL4I1_shAHR_exp_data_variables.txt",sep = "\t", stringsAsFactors = F)

setwd("../Results/IL4I1_microarrays/IL4I1_shAHR/")
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
sampleNames(raws_cel) <- s_names

#### Create phenodata & metadata ####
condition <- exp_data_variables$group
sample_id <- exp_data_variables$sample_id
date_seed <- as.factor(exp_data_variables$date_seed)
who_seed <- exp_data_variables$who_seed
who_RNA <- exp_data_variables$who_RNA
who_cDNA <- exp_data_variables$who_cDNA
who_qPCR <- exp_data_variables$who_qPCR
info <- data.frame(sample_id=sample_id, condition=condition, date_seed=date_seed, who_seed=who_seed, who_RNA=who_RNA,
                   who_cDNA=who_cDNA, who_qPCR=who_qPCR)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'date_seed', 'who_seed', 'who_RNA', 'who_cDNA', 'who_qPCR'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd

rma_2sd <- gen_rma_FUN(raws_cel, sd_val = 2, snames = s_names, ct_off = 0.85)

#### DGE ####
# Design matrix
design_eset <- model.matrix(~0+condition+who_seed, data=pd@data)
colnames(design_eset) <- gsub("condition","",colnames(design_eset))
arr_wts <- arrayWeights(oligo::exprs(rma_2sd), design = design_eset)
contrast_matrix <- makeContrasts(IL4I1_only=IL4I1_em_c-em_IL4I1_em_c,
                                 IL4I1_AHR_Dep_SH1= IL4I1_sh1-IL4I1_em_c,
                                 SH1_no_IL4I1= em_IL4I1_sh1-em_IL4I1_em_c,
                                 IL4I1_AHR_Dep_SH2=IL4I1_sh2-IL4I1_em_c,
                                 SH2_no_IL4I1=em_IL4I1_sh2-em_IL4I1_em_c,
                                 levels = design_eset)
## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset, weights = arr_wts)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_IL4I1_shAHR_u87.rds")

## all tts
tt_2sd_all <- topTable(fit_2sd_eB, number = Inf,adjust.method = "BH")
tt_2sd_u_all <- tt_2sd_all[,c(29,25,40:48)]

tts_2sd_u <- map(1:5, function(i){
  tt <- topTable(fit_2sd_eB, coef = i,  sort.by = "M",  number = Inf, adjust.method = "BH")
  tt <- tt[,c(29,25,40:45)]
})

names(tts_2sd_u) <- colnames(contrast_matrix)

## annotation update
library(dplyr)
tts_2sd_u$IL4I1_only <- annot_hgcn_FUN(tts_2sd_u$IL4I1_only, hgnc = hgnc)
tts_2sd_u$IL4I1_AHR_Dep_SH1 <- annot_hgcn_FUN(tts_2sd_u$IL4I1_AHR_Dep_SH1, hgnc = hgnc)
tts_2sd_u$SH1_no_IL4I1 <- annot_hgcn_FUN(tts_2sd_u$SH1_no_IL4I1, hgnc = hgnc)
tts_2sd_u$IL4I1_AHR_Dep_SH2 <- annot_hgcn_FUN(tts_2sd_u$IL4I1_AHR_Dep_SH2, hgnc = hgnc)
tts_2sd_u$SH2_no_IL4I1 <- annot_hgcn_FUN(tts_2sd_u$SH2_no_IL4I1, hgnc = hgnc)

## add AHR signature genes
tts_2sd_u_AHR <- map(tts_2sd_u, function(x,y){
  x$AHR_target <- match(x$hGene, y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x
}, y=AHR_genes$Gene)

saveRDS(tts_2sd_u_AHR,"./RDS/tts_2sd_u_AHR.rds")

walk2(tts_2sd_u_AHR, paste("./TopTables/tt_2sd_",names(tts_2sd_u),"_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")

#### AHR signature enrichment
roast_FUN <- function(tt_ids, cont_mat,gs,rma_obj, des_mat, arr_wts){
AHR_indexed_hi_lo <- ids2indices(gene.sets = gs, identifiers = tt_ids[,1])
prbsets_gsa_2sd <- rownames(tt_ids)
gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_obj@featureData@data$probesetid)
gsa_2sd_eset <- rma_obj[gsa_2sd_idx]
AHR_hi_lo_roast <- roast(oligo::exprs(gsa_2sd_eset),weights=arr_wts,index = AHR_indexed_hi_lo, design = des_mat,
                             contrast = cont_mat)
AHR_hi_lo_roast
}

all_roast <- map2(tts_2sd_u_AHR, as.data.frame(contrast_matrix), roast_FUN, gs=AHR_genes$Gene, rma_obj=rma_2sd, des_mat=design_eset,arr_wts=arr_wts) %>% do.call(rbind,.)
all_roast <- data.frame(Condition=rownames(all_roast), all_roast, stringsAsFactors = F)
write.table(all_roast, "./GSA/AHR_enrichment_AHR_KD_all.txt", row.names = F, sep = "\t")

## barcodeplots
pdf("./Figures/barcodeplot_m_t_stat_AHR_KD_all.pdf", width = 8, height = 12)
par(mfrow=c(3,2))
barcodeplot(tts_2sd_u_AHR$IL4I1_only$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$IL4I1_only$hGene)[[1]],main=all_roast$Condition[1])
barcodeplot(tts_2sd_u_AHR$IL4I1_AHR_Dep_SH1$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$IL4I1_AHR_Dep_SH1$hGene)[[1]],main=all_roast$Condition[2])
barcodeplot(tts_2sd_u_AHR$SH1_no_IL4I1$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$SH1_no_IL4I1$hGene)[[1]],main=all_roast$Condition[3])
barcodeplot(tts_2sd_u_AHR$IL4I1_AHR_Dep_SH2$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$IL4I1_AHR_Dep_SH2$hGene)[[1]],main=all_roast$Condition[4])
barcodeplot(tts_2sd_u_AHR$SH2_no_IL4I1$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$SH2_no_IL4I1$hGene)[[1]],main=all_roast$Condition[5])
dev.off()