#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the analysis of the U87MG cells treated with PP, HPP and I3P.
#################################################
library(AffyGEx)
library(ggrepel)

data("hgnc","msig.data.lists")

## add AHR signature genes
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", sep = "\t")

#### loading experimental condition covariates ####
exp_data_variables <- read.delim("../IL4I1_scripts_microarrays_metadata/metabolites_metadata.txt", sep = "\t", stringsAsFactors = F)

setwd("../Results/IL4I1_microarrays/Metabolites/")
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
condition <- gsub("(\\U87_)","",s_names) %>% gsub("_0.","",.) %>% gsub("DMSO","ctrl",.) %>% gsub("40uM","",.)
sample_id <- c(paste0("ctrl_",1:5), paste0("HPP_",1:5), paste0("I3P_",1:5),paste0("PP_",1:5))
date_of_exp <- as.factor(exp_data_variables$Date.performed)
who_person <- exp_data_variables$Whow.did.it.
rna_person <- exp_data_variables$RNA
info <- data.frame(sample_id=sample_id, condition=condition, date=date_of_exp, who=who_person, rna=rna_person)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'date', 'who', 'rna'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd

rma_2sd <- gen_rma_FUN(raws_cel, sd_val = 2, snames = s_names, ct_off = 0.25)

#### DGE ####
# Design matrix
design_eset <- model.matrix(~0+condition+date, data=pd@data)
colnames(design_eset) <- gsub("condition","",colnames(design_eset))
contrast_matrix <- makeContrasts(HPP=HPP-ctrl, I3P=I3P-ctrl, PP=PP-ctrl, levels = design_eset)

## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_metabolites_u87.rds")

tts_2sd_u <- map(colnames(contrast_matrix), function(i){
  tt <- topTable(fit_2sd_eB, coef = i,  sort.by = "M",  number = Inf, adjust.method = "BH")
  tt <- tt[,c(29,25,40:45)]
})

names(tts_2sd_u) <- colnames(contrast_matrix)

## annotation update
library(dplyr)
tts_2sd_u$HPP <- annot_hgcn_FUN(tts_2sd_u$HPP, hgnc = hgnc)
tts_2sd_u$I3P <- annot_hgcn_FUN(tts_2sd_u$I3P, hgnc = hgnc)
tts_2sd_u$PP <- annot_hgcn_FUN(tts_2sd_u$PP, hgnc = hgnc)

tts_2sd_u_AHR <- map(tts_2sd_u, function(x,y){
  x$AHR_target <- match(x$hGene, y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x
}, y=AHR_genes$Gene)

walk2(tts_2sd_u_AHR, paste("./TopTables/tt_2sd_",names(tts_2sd_u),"_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")

## Volcanoplots for comparison between the metabolites
tt_AHR_FUN <- function(x,y,z,a){
  x$AHR_target <- match(x[,1], y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x$logPV <- -log10(x$P.Value)
  x$colorcode <- "grey90"
  x$colorcode[which(x$logFC <=(-1*z) & x$logPV >=a)] <- "grey85"
  x$colorcode[which(x$logFC >=z & x$logPV >=a)] <- "grey85"
  x$colorcode[which(x$AHR_target == "yes" & x$logFC >=z & x$logPV >=a)] <- "red"
  x$colorcode[which(x$AHR_target == "yes" & x$logFC <=(-1*z) & x$logPV >=a)] <- "red" 
  x
}

tt_I3P <- tt_AHR_FUN(tts_2sd_u$I3P, AHR_genes$Gene, 0.58, 2)

v_p_I3P <- ggplot(tt_I3P, aes(x=logFC, y=logPV, color=colorcode, label=hGene))+xlim(c(-1.5,2))+ylim(c(0,12.5))+
  geom_point(color="grey90")+geom_point(data = tt_I3P[which(tt_I3P$AHR_target=="yes"),],size=5,show.legend = F,
                                        mapping = aes(x=logFC, y=tt_I3P$logPV[which(tt_I3P$AHR_target=="yes")]))+
  geom_label_repel(data = tt_I3P[which(tt_I3P$colorcode=="red"),],show.legend = F,
                   color="black", label.size = 0,
                   point.padding = 0.15,
                   segment.color = 'red')+
  theme_bw()+geom_vline(xintercept = 0.58, linetype="dashed", color = "grey")+
  geom_vline(xintercept = -0.58, linetype="dashed", color = "grey")+
  geom_hline(yintercept=2, linetype="dashed", color = "grey")+
  theme(legend.position = "right",
        axis.text = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),panel.grid = element_blank())+
  scale_color_manual(values = c("grey65","#E41A1C"))


tt_HPP <- tt_AHR_FUN(tts_2sd_u$HPP, AHR_genes$Gene, 0.58, 2)

v_p_HPP <- ggplot(tt_HPP, aes(x=logFC, y=logPV, color=colorcode, label=hGene))+xlim(c(-1.5,2))+ylim(c(0,12.5))+
  geom_point(color="grey90")+geom_point(data = tt_HPP[which(tt_HPP$AHR_target=="yes"),],size=5,show.legend = F,
                                        mapping = aes(x=logFC, y=tt_HPP$logPV[which(tt_HPP$AHR_target=="yes")]))+
  geom_label_repel(data = tt_HPP[which(tt_HPP$colorcode=="red"),],show.legend = F,
                   color="black", label.size = 0,
                   point.padding = 0.15,
                   segment.color = 'red')+
  theme_bw()+geom_vline(xintercept = 0.58, linetype="dashed", color = "grey")+
  geom_vline(xintercept = -0.58, linetype="dashed", color = "grey")+
  geom_hline(yintercept=2, linetype="dashed", color = "grey")+
  theme(legend.position = "right",
        axis.text = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),panel.grid = element_blank())+
  scale_color_manual(values = c("grey65","#E41A1C"))

tt_PP <- tt_AHR_FUN(tts_2sd_u$PP, AHR_genes$Gene, 0.58, 2)

v_p_PP <- ggplot(tt_PP, aes(x=logFC, y=logPV, color=colorcode, label=hGene))+xlim(c(-1.5,2))+ylim(c(0,12.5))+
  geom_point(color="grey90")+geom_point(data = tt_PP[which(tt_PP$AHR_target=="yes"),],size=5,show.legend = F,
                                        mapping = aes(x=logFC, y=tt_PP$logPV[which(tt_PP$AHR_target=="yes")]))+
  geom_label_repel(data = tt_PP[which(tt_PP$colorcode=="red"),],show.legend = F,
                   color="black",label.size = 0,
                   point.padding = 0.15,
                   segment.color = 'red')+
  theme_bw()+geom_vline(xintercept = 0.58, linetype="dashed", color = "grey")+
  geom_vline(xintercept = -0.58, linetype="dashed", color = "grey")+
  geom_hline(yintercept=2, linetype="dashed", color = "grey")+
  theme(legend.position = "right",
        axis.text = element_text(size = 18, color = "black"),
        text = element_text(size = 18, color = "black"),panel.grid = element_blank())+
  scale_color_manual(values = c("grey65","#E41A1C"))

pdf("./Figures/I3P_volcano_plot.pdf", width = 8, height = 8, pointsize = 28)
v_p_I3P
dev.off()

pdf("./Figures/HPP_volcano_plot.pdf", width = 8, height = 8, pointsize = 28)
v_p_HPP
dev.off()

pdf("./Figures/PP_volcano_plot.pdf", width = 8, height = 8, pointsize = 28)
v_p_PP
dev.off()
