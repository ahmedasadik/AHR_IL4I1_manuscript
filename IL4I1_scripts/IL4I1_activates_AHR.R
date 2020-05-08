#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes how we came across IL4I1 as an enzyme producing ligands activating AHR.
#################################################

## Load libraries
library(purrr)
library(fmsb)

## Source the functions and parameters files
source("../functions_and_parameters.R")

# Selected enzymes of the Tryptophan degradation pathway
trp_enz_sel <- c("IL4I1","IDO1","TDO2","TPH1","TPH2","DDC", "IDO2")

# TCGA_voom
TCGA_voom <- readRDS("../Zenodo_download/TCGA_DGE_voom_annot.rds")

# TCGA tumors and the modules correlating with AHR (this shows both positive and negative correlations - pearson)
TCGA_MEs_GSVA_cors <- readRDS("../Results/RDS/TCGA_MEs_GSVA_cors.rds")

# Modules correlating only positively with AHR activation
modules_with_enzymes_pos <- readRDS("../Results/RDS/modules_with_enzymes_pos.rds")

# TCGA_GSVA
TCGA_gsva <- readRDS("../Zenodo_download/TCGA_GSVA_scores_safely.rds")
TCGA_GSVA <- TCGA_gsva[tcga_names]
TCGA_GSVA <- map(TCGA_GSVA, function(x){x$result})

# TCGA_modules_genes
TCGA_modules_genes <- readRDS("../Zenodo_download/TCGA_modules_genes.rds")

# TCGA_MEs_GSVA_cors_pos_only_overlap_GT
TCGA_MEs_GSVA_cors_pos_only_overlap_GT <- readRDS("../Results/RDS/TCGA_MEs_GSVA_cors_pos_only_overlap_GT.rds")

#######################################################################
## Since IDO1 and TDO2 are not driving AHR activation in all tumors, ##
## what is the status of other enzymes of Trp degradation pathway?   ##
#######################################################################
## The incidence plots of selected enzymes from Trp degradation pathway
final_TCGA_MEs_GSVA_corrs_AHR_enzymes <- map2(TCGA_MEs_GSVA_cors, modules_with_enzymes_pos, function(cors, enzs){
  cors$enz_aaas_l <- ""
  cors$enz_ahr_l <- ""
  cors$enz_aaas_g <- ""
  cors$enz_ahr_g <- ""
  cors$enz_aaas_ahr <- ""
  enz_modules <- enzs$module
  for (i in 1:length(cors$modules)){
    enz_mod_idx <- which(enz_modules==cors$modules[i])
    aaas_idx <- which(enzs[enz_mod_idx,"aaa_enzymes"]=="yes")
    ahr_idx <- which(enzs[enz_mod_idx,"AHR_genes"]=="yes")
    cors$enz_aaas_l[i] <- length(aaas_idx)
    cors$enz_ahr_l[i] <- length(ahr_idx)
    cors$enz_aaas_g[i] <- paste(enzs$genes[enz_mod_idx][aaas_idx],collapse = ";")
    cors$enz_ahr_g[i] <- paste(enzs$genes[enz_mod_idx][ahr_idx],collapse = ";")
    cors$enz_aaas_ahr[i] <- paste(enzs$genes[enz_mod_idx][intersect(aaas_idx, ahr_idx)],collapse = ";")
  }
  cors
})

names(final_TCGA_MEs_GSVA_corrs_AHR_enzymes) <- tcga_names

final_TCGA_MEs_GSVA_corrs_AHR_enzymes <- map2(final_TCGA_MEs_GSVA_corrs_AHR_enzymes, tcga_names, function(x,y){
  df <- data.frame(tumor= gsub("TCGA_","",y),x)
}) %>% do.call(rbind,.)

final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IL4I1 <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IDO1 <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$TDO2 <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$none <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IL4I1[grep("IL4I1",final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_g)] <- "yes"
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IDO1[grep("IDO1",final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_g)] <- "yes"
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$TDO2[grep("TDO2",final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_g)] <- "yes"
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$none[which(final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_l==0 & final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_ahr_l==0)] <- "ns"

## dataframe showing the incidence for each individual tumor
enz_no_in_modules_per_tumor <- map(gsub("TCGA_","",tcga_names), function(a,b,c){
  tumor_df <- b[b$tumor==a,]
  incid_df <- map2(c,1, function(x,y){
    if(length(grep(x,tumor_df$enz_aaas_g))==0){
      0
    } else {
      y
    }
    }) %>% do.call("c",.)
    },
  b=final_TCGA_MEs_GSVA_corrs_AHR_enzymes, c=trp_enz_sel)

names(enz_no_in_modules_per_tumor) <- gsub("TCGA_","",tcga_names)
enz_no_in_modules_per_tumor_df <- do.call(rbind,enz_no_in_modules_per_tumor)
colnames(enz_no_in_modules_per_tumor_df) <- trp_enz_sel

write.table(enz_no_in_modules_per_tumor_df, "../Results/Tables/enz_no_in_modules_per_tumor_df.txt", sep = "\t")

radar_enz_no_in_modules_per_tumor_df <- rbind(rep(1,length(tcga_names)),rep(0,length(tcga_names)),t(enz_no_in_modules_per_tumor_df))
rownames(radar_enz_no_in_modules_per_tumor_df)[1:2] <- c("low","high")
new_order <- c("low", "high", "IL4I1", "IDO1", "IDO2", "TDO2", "DDC", "TPH1", "TPH2")
radar_enz_no_in_modules_per_tumor_df <- radar_enz_no_in_modules_per_tumor_df[new_order,]

write.table(radar_enz_no_in_modules_per_tumor_df, "../Results/Tables/radar_enz_no_in_modules_per_tumor_df.txt", sep = "\t")

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9),
                 rgb(0.3,0.6,0.5,0.9), rgb(0.4,0.7,0.3,0.9) , rgb(0.6,0.5,0.8,0.9))
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4),
             rgb(0.3,0.6,0.5,0.4), rgb(0.4,0.7,0.3,0.4) , rgb(0.6,0.5,0.8,0.4))

pdf("../Results/Figures/TRP_enz_ME_incidences_pos_cor_per_tumor_clockwise2.pdf", height = 12, width = 16, pointsize = 16)
map2(new_order[-c(1:2)],c(2,rep(1,6)), function(x,i){
  radarchart(as.data.frame(radar_enz_no_in_modules_per_tumor_df[c("low","high",x),]), axistype=0 ,title = x,
             #custom polygon
             pcol=colors_border[i] , pfcol=colors_in[i] , plwd=4 , plty=1,
             #custom the grid
             cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
             #custom labels
             vlcex=0.8 
  )
})
dev.off()

########################################################
## The AHR biological functions associated with IL4I1 ##
########################################################

## The overlapping modules between the cor and GT
pos_only_mods <- map2(TCGA_MEs_GSVA_cors_pos_only_overlap_GT,tcga_names, function(x,y){
  paste(gsub("TCGA_","",y),x$modules,sep = "_")
  }) %>%  unlist(.)

## Read the wgcna_ontology files
wgcna_ontologies <- read.csv("../Resources/gene_ontologie_matrix_all_modules_WGCNA_AHR_functions.csv",stringsAsFactors = F)
wgcna_ontologies$comb <-paste(wgcna_ontologies$Tumor,wgcna_ontologies$modules,sep = "_") %>% gsub(" ","",.)

wgcna_ontologies <- wgcna_ontologies[wgcna_ontologies$comb%in%pos_only_mods,]

## The Gene ontology bargraph
ontologies <- c("Angiogenesis",	"Drug metabolism",	"External stress response",	"Hemopoiesis",	"Immune modulation",	"Lipid metabolsim",	"Motility")

## The modules that have IL4I1 in all tumors
IL4I1_modules <- map(TCGA_modules_genes,function(x){x$module[x$genes=="IL4I1"]})

IL4I1_modules2 <- map(IL4I1_modules, function(x){
  if(length(x)>0){
    x
  } else {
    "empty"
  }
})
IL4I1_modules2 <- paste(gsub("TCGA_","",tcga_names), unlist(IL4I1_modules2),sep = "_")
ont_incid_IL4I1 <- wgcna_ontologies[wgcna_ontologies$comb%in%IL4I1_modules2,]%>% .[,3:9] %>% sapply(.,as.numeric)%>%colSums()
ont_df_IL4I1 <- data.frame(ontologies,ont_incid_IL4I1,stringsAsFactors = F)

# Create data: note in High school for Jonathan:
radar_df_IL4I1 <- data.frame(rep(max(ont_df_IL4I1$ont_incid_IL4I1), length(ont_df_IL4I1$ont_incid_IL4I1)),
                             rep(0, length(ont_df_IL4I1$ont_incid_IL4I1)),ont_df_IL4I1$ont_incid_IL4I1)%>%t()
colnames(radar_df_IL4I1) <- ont_df_IL4I1$ontologies

pdf("../Results/Figures/Incidince_AHR_ontology groups_associated_with_IL4I1.pdf",width = 12,height = 8)
radarchart(as.data.frame(radar_df_IL4I1[,order(radar_df_IL4I1[3,])]), axistype=0 ,title = "IL4I1 ontology groups",
           #custom polygon
           pcol=colors_border[2] , pfcol=colors_in[2] , plwd=4 , plty=1,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
           #custom labels
           vlcex=0.8 
)
dev.off()