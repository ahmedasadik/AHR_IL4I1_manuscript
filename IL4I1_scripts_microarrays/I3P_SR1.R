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
exp_data_variables <- read.delim("../IL4I1_scripts_microarrays_metadata/I3P_SR1_exp_data_variables.txt", sep = "\t", stringsAsFactors = F)

setwd("../Results/IL4I1_microarrays/I3P_SR1/")
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
s_names <- gsub("\\_\\HuGene\\-2\\_0\\-st\\_\\.CEL","",sampleNames(raws_cel)) %>%  gsub("GSM.([0-9]*_)","",.) %>% gsub("-_","",.)
sampleNames(raws_cel) <- s_names

#### Create phenodata & metadata ####
condition <- c(rep(c("ctrl", "I3P", "I3P_SR1", "SR1"),2), rep(c("ctrl", "I3P", "SR1", "I3P_SR1"), 2))
sample_id <- exp_data_variables$Sample
date_exp <- rep(c("d1", "d2", "d3", "d4"), each=4)
conc_rna <- exp_data_variables$Conc..µg.µL.
info <- data.frame(sample_id=sample_id, condition=condition, date_exp=date_exp, conc=conc_rna)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'date_exp', 'conc'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd

rma_2sd <- gen_rma_FUN(raws_cel, sd_val = 2,snames = s_names,ct_off = 0.9)

#### DGE ####
# Design matrix
design_eset <- model.matrix(~0+condition+date_exp, data=pd@data)
colnames(design_eset) <- gsub("condition","",colnames(design_eset))
design_eset <- cbind(design_eset, as.character(pd@data$conc)%>%as.numeric(.))
colnames(design_eset)[8] <- "conc"
arr_wts <- arrayWeights(oligo::exprs(rma_2sd), design = design_eset)
contrast_matrix <- makeContrasts(I3P=I3P-ctrl, SR1=SR1-ctrl, I3P_SR1=I3P_SR1-ctrl, Diff=(I3P-ctrl)-(I3P_SR1-ctrl),
                                 levels = design_eset)
## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset, weights = arr_wts)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_migration.rds")

## all tts
tt_2sd_all <- topTable(fit_2sd_eB, number = Inf,adjust.method = "BH")
tt_2sd_u_all <- tt_2sd_all[,c(29,25,40:45)]

tts_2sd_u <- map(1:4, function(i){
  tt <- topTable(fit_2sd_eB, coef = i,  sort.by = "M",  number = Inf, adjust.method = "BH")
  tt <- tt[,c(29,25,40:45)]
})

names(tts_2sd_u) <- colnames(contrast_matrix)

## annotation update
library(dplyr)
tts_2sd_u[[1]] <- annot_hgcn_FUN(tts_2sd_u[[1]], hgnc = hgnc)
tts_2sd_u[[2]] <- annot_hgcn_FUN(tts_2sd_u[[2]], hgnc = hgnc)
tts_2sd_u[[3]] <- annot_hgcn_FUN(tts_2sd_u[[3]], hgnc = hgnc)
tts_2sd_u[[4]] <- annot_hgcn_FUN(tts_2sd_u[[4]], hgnc = hgnc)

## add AHR signature genes
tts_2sd_u_AHR <- map(tts_2sd_u, function(x,y){
  x$AHR_target <- match(x$hGene, y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x
}, y=AHR_genes$Gene)

walk2(tts_2sd_u_AHR, paste("./TopTables/tt_2sd_",names(tts_2sd_u),"_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")

#### GSA ####
gsa_res <- gsa_FUN(top_t = tts_2sd_u_AHR$I3P, rma_obj = rma_2sd, exprmnt = "migration", cntrst = "many",
                   cnt_mat = contrast_matrix, res_path = "./GSA/",wts = arr_wts,
                   msig.data.lists = msig.data.lists, msigs = 14, pltfrm = "Affy", d_m = design_eset)

names(gsa_res) <- names(tts_2sd_u_AHR)

for(i in 1:length(gsa_res)){
  walk(1:length(gsa_res[[i]]), safely(gsa_bar_plot_FUN), gsa_ls=gsa_res[[i]], res_path="/Figures/",
       coi=names(gsa_res)[i], wd=12, ht=8, file_type="pdf")
}

gsa_paths <- list.files("./GSA/",pattern = "_migration_",full.names = T)
gsa_paths <- gsa_paths[grep("\\.txt",gsa_paths)] %>% .[-1]
All_tts_GSAs <- map(gsa_paths,read.delim,stringsAsFactors=F)
gsa_paths_names <- gsub("\\.//|\\.txt","",gsa_paths)
names(All_tts_GSAs) <- gsa_paths_names

condition_names<-gsub("gsa_|_migration.*.","",gsa_paths_names)

All_tts_GSAs_05 <- map2(All_tts_GSAs,condition_names, function(x,y){
  df <- data.frame(condition=y, pathways=rownames(x),x,stringsAsFactors = F)
  df <- df[df$PValue<=0.05,]
  df$Direction <- factor(df$Direction, levels = c("Up","Down"))
  df$logPV <- -log10(df$PValue)
  if(length(df$pathways)>30){
    df[1:15,]
  } else {
    df
  }
})

names(All_tts_GSAs_05) <- gsub("\\./GSA//","",names(All_tts_GSAs_05))

## Heatmap like structure
require(corrplot)
# Hallmarks
hallmarks_df <- dplyr::full_join(All_tts_GSAs_05$gsa_I3P_migration_res_hallmark,
                                 All_tts_GSAs_05$gsa_I3P_SR1_migration_res_hallmark,
                                 by=c("pathways"))
hallmarks_df_r <- data.frame(I3P=hallmarks_df$Direction.x, I3P_SR1=hallmarks_df$Direction.y,stringsAsFactors = F)
hallmarks_df_r$I3P <- ifelse(hallmarks_df_r$I3P=="Up",1,ifelse(hallmarks_df_r$I3P=="Down",-1,0))
hallmarks_df_r$I3P_SR1 <- ifelse(hallmarks_df_r$I3P_SR1=="Up",1,ifelse(hallmarks_df_r$I3P_SR1=="Down",-1,0))

hallmarks_df_p <- data.frame(I3P=hallmarks_df$NGenes.x, I3P_SR1=hallmarks_df$NGenes.y,stringsAsFactors = F)
hallmarks_df_p$I3P <- ifelse(is.na(hallmarks_df_r$I3P),0,hallmarks_df_p$I3P) 
hallmarks_df_p$I3P <- ifelse(hallmarks_df_r$I3P>0, hallmarks_df_p$I3P,ifelse(hallmarks_df_r$I3P<0,(hallmarks_df_p$I3P *(-1)),hallmarks_df_p$I3P))
hallmarks_df_p$I3P_SR1 <- ifelse(is.na(hallmarks_df_r$I3P_SR1),0,hallmarks_df_p$I3P_SR1) 
hallmarks_df_p$I3P_SR1 <- ifelse(hallmarks_df_r$I3P_SR1>0, hallmarks_df_p$I3P_SR1,ifelse(hallmarks_df_r$I3P_SR1<0,(hallmarks_df_p$I3P_SR1 *(-1)),hallmarks_df_p$I3P_SR1))

rownames(hallmarks_df_r) <- rownames(hallmarks_df_p) <- hallmarks_df$pathways

pdf("./Figures/hallmarks_heatmap.pdf", height = 8,width = 12)
corrplot(corr = as.matrix(hallmarks_df_r),p.mat = as.matrix(hallmarks_df_r),method = "square",tl.col = "black",cl.pos = "n",
         is.corr = FALSE,insig = "n",na.label.col = "white",col = c("blue","red"),mar = c(0,0,0.5,12),
)
dev.off()

# KEGG
kegg_df <- dplyr::full_join(All_tts_GSAs_05$gsa_I3P_migration_res_kegg,
                            All_tts_GSAs_05$gsa_I3P_SR1_migration_res_kegg,
                            by=c("pathways"))

kegg_df_r <- data.frame(I3P=kegg_df$Direction.x, I3P_SR1=kegg_df$Direction.y,stringsAsFactors = F)
kegg_df_r$I3P <- ifelse(kegg_df_r$I3P=="Up",1,ifelse(kegg_df_r$I3P=="Down",-1,0))
kegg_df_r$I3P_SR1 <- ifelse(kegg_df_r$I3P_SR1=="Up",1,ifelse(kegg_df_r$I3P_SR1=="Down",-1,0))

kegg_df_p <- data.frame(I3P=kegg_df$NGenes.x, I3P_SR1=kegg_df$NGenes.y,stringsAsFactors = F)
kegg_df_p$I3P <- ifelse(is.na(kegg_df_r$I3P),0,kegg_df_p$I3P) 
kegg_df_p$I3P <- ifelse(kegg_df_r$I3P>0, kegg_df_p$I3P,ifelse(kegg_df_r$I3P<0,(kegg_df_p$I3P *(-1)),kegg_df_p$I3P))
kegg_df_p$I3P_SR1 <- ifelse(is.na(kegg_df_r$I3P_SR1),0,kegg_df_p$I3P_SR1) 
kegg_df_p$I3P_SR1 <- ifelse(kegg_df_r$I3P_SR1>0, kegg_df_p$I3P_SR1,ifelse(kegg_df_r$I3P_SR1<0,(kegg_df_p$I3P_SR1 *(-1)),kegg_df_p$I3P_SR1))

rownames(kegg_df_r) <- rownames(kegg_df_p) <- kegg_df$pathways

pdf("./Figures/kegg_heatmap.pdf", height = 8,width = 12)
corrplot(corr = as.matrix(kegg_df_r),p.mat = as.matrix(kegg_df_r),method = "square",tl.col = "black",cl.pos = "n",
         is.corr = FALSE,insig = "n",na.label.col = "white",col = c("blue","red"),mar = c(0,0,0.5,12),tl.cex = 0.8,
)
dev.off()


#### AHR signature enrichment
roast_FUN <- function(tt_ids, cont_mat,gs,rma_obj, des_mat, arr_wts){
  AHR_indexed_hi_lo <- ids2indices(gene.sets = gs, identifiers = tt_ids[,1])
  prbsets_gsa_2sd <- rownames(tt_ids)
  gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_obj@featureData@data$probesetid)
  gsa_2sd_eset <- rma_obj[gsa_2sd_idx]
  AHR_hi_lo_roast <- roast(oligo::exprs(gsa_2sd_eset),weights=arr_wts,index = AHR_indexed_hi_lo, design = des_mat,
                               contrast = cont_mat, set.statistic = "floormean")
  AHR_hi_lo_roast
}

all_roast <- map2(tts_2sd_u_AHR, as.data.frame(contrast_matrix), roast_FUN, gs=AHR_genes$Gene, rma_obj=rma_2sd, des_mat=design_eset, arr_wts=arr_wts) %>% do.call(rbind,.)
all_roast <- data.frame(Condition=rownames(all_roast), all_roast, stringsAsFactors = F)

write.table(all_roast, "./GSA/AHR_enrichment_migration_all.txt", row.names = F, sep = "\t")

pdf("./Figures/barcodeplot_I3P_SR1.pdf", width = 12, height = 8)
barcodeplot(tts_2sd_u_AHR[[1]]$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR[[1]]$hGene)[[1]],main=all_roast$Condition[1])
barcodeplot(tts_2sd_u_AHR[[2]]$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR[[2]]$hGene)[[1]],main=all_roast$Condition[2])
barcodeplot(tts_2sd_u_AHR[[3]]$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR[[3]]$hGene)[[1]],main=all_roast$Condition[3])
barcodeplot(tts_2sd_u_AHR[[4]]$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR[[4]]$hGene)[[1]],main=all_roast$Condition[4])
dev.off()
