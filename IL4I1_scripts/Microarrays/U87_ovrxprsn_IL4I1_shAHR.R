# This is script is for the analysis of the microarray of U87 IL4I1 overexpressing cells with or without shAHR (1 and 2)
library(AffyPipelineHG2ST)

data("TFSBDB_list","TFSBDB_list_comb","TFSBDB_list_comb_genes",
     "TFSBDB_list_genes","hgnc","msig.data.lists")

setwd("./IL4I1/Results_U87_IL4I1_ovrxprsn_AHR_KD/")

dir.create("./RDS")
dir.create("./TopTables")
dir.create("./GSA")
dir.create("./Figures")
dir.create("./TFSDB")

#### generate raw_cel dataset ####
raw_path <- "/home/data/Raws/IL4I1/shAHR/"
list_cels <- list.files(raw_path, pattern = "CEL$", full.names = T)
raws_cel <- read.celfiles(list_cels)
saveRDS(raws_cel, "./RDS/raws_cel.rds")

#### Assigning sample names ####
s_names <- gsub("\\_\\(HuGene\\-2\\_0\\-st\\)\\.CEL","",sampleNames(raws_cel))
sampleNames(raws_cel) <- s_names

#### loading experimental condition covariates ####
exp_data_variables <- xlsx::read.xlsx("/home/data/Raws/IL4I1/20180206 U87-IL4I1-shAHR_RNA microarray-1.xlsx",1, stringsAsFactors=F)
exp_data_variables <- exp_data_variables[,-13]
colnames(exp_data_variables) <- c("group", "sample_name", "sample_id", "conc", "RNA", "H2O", "date_seed", "who_seed", "passage",
                                  "who_RNA","who_cDNA","who_qPCR")
exp_data_variables <- exp_data_variables[-c(1,2,21),]

exp_data_variables <- exp_data_variables[c(2,8,14,3,9,15,1,7,13,5,11,17,6,12,18,4,10,16),]
exp_data_variables$group <- gsub("_em_sh1","_sh1",exp_data_variables$group) %>% gsub("_em_sh2","_sh2",.)

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

rma_2sd <- gen_rma_cel_nocross_median_FUN(raws_cel, sd_val = 2,snames = s_names,ct_off = 0.85)

### PCA plots to determine possible batch effects raws ####
pca_raw <- prcomp(t(exprs(rma(raws_cel))), scale. = T)
pdf("./Figures/pca_raws_IL4I1_sh.pdf", width = 6, height = 6)
factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA IL4I1 shAHR groups", habillage=pd@data$condition, pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$date_seed))],
                 title = "PCA IL4I1 shAHR date seed", habillage=exp_data_variables$date_seed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$passage))],
                 title = "PCA IL4I1 shAHR passage", habillage=exp_data_variables$passage,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$who_seed))],
                 title = "PCA IL4I1 shAHR who seed", habillage=exp_data_variables$who_seed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$who_RNA))],
                 title = "PCA IL4I1 shAHR who RNA", habillage=exp_data_variables$who_RNA,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$who_cDNA))],
                 title = "PCA IL4I1 shAHR who cDNA", habillage=exp_data_variables$who_cDNA,pointsize = 6, mean.point = FALSE)
dev.off()

### PCA plots to determine possible batch effects ####
pca_rma <- prcomp(t(exprs(rma_2sd)), scale. = T)
pdf("./Figures/pca_normalized_IL4I1 shAHR.pdf", width = 6, height = 6)
factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA IL4I1 shAHR groups", habillage=pd@data$condition, pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$date_seed))],
                 title = "PCA IL4I1 shAHR date seed", habillage=exp_data_variables$date_seed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$passage))],
                 title = "PCA IL4I1 shAHR passage", habillage=exp_data_variables$passage,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$who_seed))],
                 title = "PCA IL4I1 shAHR who seed", habillage=exp_data_variables$who_seed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$who_RNA))],
                 title = "PCA IL4I1 shAHR who RNA", habillage=exp_data_variables$who_RNA,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$who_cDNA))],
                 title = "PCA IL4I1 shAHR who cDNA", habillage=exp_data_variables$who_cDNA,pointsize = 6, mean.point = FALSE)
dev.off()

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
walk2(tts_2sd_u, paste("./TopTables/tt_2sd_",names(tts_2sd_u),".txt",sep = ""), write.table, quote = F, sep = "\t")

## add AHR signature genes
AHR_genes <- read.delim("/home/analyses/Projects_analyses/AHR_IL4I1/AHR/Results/Signature/overlapping_AHR_signature.txt", sep = "\t")
tts_2sd_u_AHR <- map(tts_2sd_u, function(x,y){
  x$AHR_target <- match(x$hGene, y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x
}, y=AHR_genes$Gene)

walk2(tts_2sd_u_AHR, paste("./TopTables/tt_2sd_",names(tts_2sd_u),"_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")

#### Combination comparisons ####
tt_all_hgene_ord <- map(tts_2sd_u, function(x){x[order(x$hGene),]})

fc_na_compilation_FUN <- function(tt, fc_ct, pv_ct){
  fc_na_FUN <- function(x, y, z, v){
    if(x >= z & y <= v){
      fc <- x
    } else if(x <= -z & y <= v){
      fc <- x
    } else {
      fc <- NA
    }
    fc
  }
  fc_values <- map2(tt$logFC, tt$P.Value, fc_na_FUN, z=fc_ct, v=pv_ct) %>% unlist()
}
fc_na_compilation_apv_FUN <- function(tt, fc_ct, pv_ct){
  fc_na_FUN <- function(x, y, z, v){
    if(x >= z & y <= v){
      fc <- x
    } else if(x <= -z & y <= v){
      fc <- x
    } else {
      fc <- NA
    }
    fc
  }
  fc_values <- map2(tt$logFC, tt$adj.P.Val, fc_na_FUN, z=fc_ct, v=pv_ct) %>% unlist()
}
fc_na_compilation_grp_FUN <-function(x){
  hi_idcs <- which(x > 0)
  lo_idcs <- which(x < 0)
  na_idcs <- which(is.na(x))
  group <- vector("character", length = length(x))
  group[hi_idcs] <- "up"
  group[lo_idcs] <- "down"
  group[na_idcs] <- "noDE"
  group
}
comp_group_assign_FUN <- function(E2, INFG, E2INFG){
  group <- c()
  if(E2 == "noDE" & INFG=="noDE" & E2INFG=="noDE"){
    group <- "noDE"
  } 
  else if(E2 != "noDE" & INFG=="noDE" & E2INFG=="noDE"){
    group <- "E2_only"
  } 
  else if (E2 == "noDE" & INFG!="noDE" & E2INFG=="noDE"){
    group <- "INFG_only"
  }
  else if (E2 != "noDE" & INFG!="noDE" & E2INFG=="noDE"){
    group <- "E2_INFG_only"
  } 
  else if (E2 == "noDE" & INFG=="noDE" & E2INFG!="noDE"){
    group <- "E2INFG_only"
  } 
  else if (E2 != "noDE" & INFG=="noDE" & E2INFG!="noDE"){
    group <- "E2_E2INFG_only"
  } 
  else if (E2 == "noDE" & INFG!="noDE" & E2INFG!="noDE"){
    group <- "INFG_E2INFG_only"
  } 
  else if (E2 != "noDE" & INFG!="noDE" & E2INFG!="noDE"){
    if(E2 == "up" & INFG=="up" & E2INFG=="up"){
      group<- "syn_up_both_up"
    } else if(E2 == "up" & INFG=="up" & E2INFG=="down"){
      group<- "syn_dn_opp_up"
    } else if(E2 == "up" & INFG=="down" & E2INFG=="up"){
      group<- "e2_up_gr_infg"
    } else if(E2 == "up" & INFG=="down" & E2INFG=="down"){
      group<- "infg_dn_gr_e2"
    } else if(E2 == "down" & INFG=="up" & E2INFG=="up"){
      group<- "infg_up_gr_e2"
    } else if(E2 == "down" & INFG=="up" & E2INFG=="down"){
      group<- "e2_dn_gr_infg"
    } else if(E2 == "down" & INFG=="down" & E2INFG=="up"){
      group<- "syn_up_opp_dn"
    } else if(E2 == "down" & INFG=="down" & E2INFG=="down"){
      group<- "syn_dn_both_dn"
    }
  } else {
    group <- NA
  }
  group
}

## PV filtering ####
FCs_comp_pv <- sapply(tt_all_hgene_ord, fc_na_compilation_FUN, fc_ct=0.3, pv_ct=0.05) %>% 
  data.frame(gene=tt_all_hgene_ord$IL4I1_only$hGene, ., stringsAsFactors=FALSE)
FCs_comp_pv <-  data.frame(FCs_comp_pv,sapply(FCs_comp_pv[,-1], fc_na_compilation_grp_FUN), stringsAsFactors = FALSE)

SH1_comb_diff_pv <- unlist(map2(FCs_comp_pv$IL4I1_only, FCs_comp_pv$IL4I1_AHR_Dep_SH1,function(x, y){
  if(is.na(x) & !is.na(y) & y > 0){
    "Induced in combination _ check IL4I1 and KD contributions"
  } else if(is.na(x) & !is.na(y) & y < 0){
    "Down regulated in combination _ check IL4I1 and KD contributions"
  } else if (!is.na(x) & is.na(y) & x>0){
    "Induced by IL4I1  and not differentially regulated in combination"
  } else if(!is.na(x) & is.na(y) & x<0){
    "Down regulated by IL4I1 and not differentially regulated in combination"
  } else if (is.na(x) & is.na(y)){
    "Check KD status"
  } else if(x > 0 & x > y){
    "Induced by IL4I1 and down regulated by KD in combination"
  } else if(x > 0 & x < y){
    "Induced by IL4I1 and even more by KD in combination"
  } else if (x < 0 & x < y){
    "Down regulated by IL4I1 and up regulated by KD in combination"
  } else if(x < 0 & x > y){
    "Down regulated by IL4I1 and even more by KD in combination"
  }
}))
SH2_comb_diff_pv <- unlist(map2(FCs_comp_pv$IL4I1_only, FCs_comp_pv$IL4I1_AHR_Dep_SH2,function(x, y){
  if(is.na(x) & !is.na(y) & y > 0){
    "Induced in combination _ check IL4I1 and KD contributions"
  } else if(is.na(x) & !is.na(y) & y < 0){
    "Down regulated in combination _ check IL4I1 and KD contributions"
  } else if (!is.na(x) & is.na(y) & x>0){
    "Induced by IL4I1  and not differentially regulated in combination"
  } else if(!is.na(x) & is.na(y) & x<0){
    "Down regulated by IL4I1 and not differentially regulated in combination"
  } else if (is.na(x) & is.na(y)){
    "Check KD status"
  } else if(x > 0 & x > y){
    "Induced by IL4I1 and down regulated by KD in combination"
  } else if(x > 0 & x < y){
    "Induced by IL4I1 and even more by KD in combination"
  } else if (x < 0 & x < y){
    "Down regulated by IL4I1 and up regulated by KD in combination"
  } else if(x < 0 & x > y){
    "Down regulated by IL4I1 and even more by KD in combination"
  }
}))

FCs_pv_df <- data.frame(FCs_comp_pv[,1:6], SH1_IL4I1_comp_group=SH1_comb_diff_pv,
                        SH2_IL4I1_comp_group=SH2_comb_diff_pv,
                        consensus_group = ifelse(SH1_comb_diff_pv == SH2_comb_diff_pv, "yes", "no"))

FCs_pv_df <- FCs_pv_df[-which(is.na(FCs_pv_df$IL4I1_only) & is.na(FCs_pv_df$IL4I1_AHR_Dep_SH1) & is.na(FCs_pv_df$SH1_no_IL4I1)),]
FCs_pv_df$IL4I1_Abs <- ifelse(!is.na(FCs_pv_df$IL4I1_only),ifelse(FCs_pv_df$IL4I1_only > 0, 2^FCs_pv_df$IL4I1_only, -1*(2^abs(FCs_pv_df$IL4I1_only))),NA)
FCs_pv_df$SH1_Abs <- ifelse(!is.na(FCs_pv_df$IL4I1_AHR_Dep_SH1),ifelse(FCs_pv_df$IL4I1_AHR_Dep_SH1 > 0, 2^FCs_pv_df$IL4I1_AHR_Dep_SH1, -1*(2^abs(FCs_pv_df$IL4I1_AHR_Dep_SH1))),NA)
FCs_pv_df$SH2_Abs <- ifelse(!is.na(FCs_pv_df$IL4I1_AHR_Dep_SH2),ifelse(FCs_pv_df$IL4I1_AHR_Dep_SH2 > 0, 2^FCs_pv_df$IL4I1_AHR_Dep_SH2, -1*(2^abs(FCs_pv_df$IL4I1_AHR_Dep_SH2))),NA)

FCs_pv_df$IL4I1_SH1_Abs_diff <- ifelse(!is.na(FCs_pv_df$IL4I1_Abs) & !is.na(FCs_pv_df$SH1_Abs),
                                       FCs_pv_df$IL4I1_Abs-FCs_pv_df$SH1_Abs, 
                                       ifelse(!is.na(FCs_pv_df$IL4I1_Abs) & is.na(FCs_pv_df$SH1_Abs),
                                              FCs_pv_df$IL4I1_Abs, 
                                              ifelse(is.na(FCs_pv_df$IL4I1_Abs) & !is.na(FCs_pv_df$SH1_Abs),
                                                     FCs_pv_df$SH1_Abs, NA)))
FCs_pv_df$IL4I1_SH2_Abs_diff <- ifelse(!is.na(FCs_pv_df$IL4I1_Abs) & !is.na(FCs_pv_df$SH2_Abs),
                                       FCs_pv_df$IL4I1_Abs-FCs_pv_df$SH2_Abs, 
                                       ifelse(!is.na(FCs_pv_df$IL4I1_Abs) & is.na(FCs_pv_df$SH2_Abs),
                                              FCs_pv_df$IL4I1_Abs, 
                                              ifelse(is.na(FCs_pv_df$IL4I1_Abs) & !is.na(FCs_pv_df$SH2_Abs),
                                                     FCs_pv_df$SH2_Abs, NA)))

FCs_pv_df$consensus_ABS <- ifelse(abs(FCs_pv_df$IL4I1_SH1_Abs_diff-FCs_pv_df$IL4I1_SH2_Abs_diff) <=0.5, "yes", "no")

FCs_pv_df$TFs <- names(TFSBDB_list_comb_genes)[match(FCs_pv_df$gene, names(TFSBDB_list_comb_genes))]
FCs_pv_df$AHR <- AHR_genes$Gene[match(FCs_pv_df$gene, AHR_genes$Gene)]
write.table(FCs_pv_df, "./TopTables/comparison_all_genes_grouping_pv.txt", row.names = F, quote = F, sep = "\t")

FCs_pv_df_AHR_only <- FCs_pv_df[!is.na(FCs_pv_df$AHR),]
write.table(FCs_pv_df_AHR_only, "./TopTables/comparison_all_genes_grouping_pv_AHR_only.txt", row.names = F, quote = F, sep = "\t")

## APV filtering ####
FCs_comp_apv <- sapply(tt_all_hgene_ord, fc_na_compilation_apv_FUN, fc_ct=0.3, pv_ct=0.05) %>% 
  data.frame(gene=tt_all_hgene_ord$IL4I1_only$hGene, ., stringsAsFactors=FALSE)
FCs_comp_apv <-  data.frame(FCs_comp_apv,sapply(FCs_comp_apv[,-1], fc_na_compilation_grp_FUN), stringsAsFactors = FALSE)

SH1_comb_diff_apv <- unlist(map2(FCs_comp_apv$IL4I1_only, FCs_comp_apv$IL4I1_AHR_Dep_SH1,function(x, y){
  if(is.na(x) & !is.na(y) & y > 0){
    "Induced in combination _ check IL4I1 and KD contributions"
  } else if(is.na(x) & !is.na(y) & y < 0){
    "Down regulated in combination _ check IL4I1 and KD contributions"
  } else if (!is.na(x) & is.na(y) & x>0){
    "Induced by IL4I1  and not differentially regulated in combination"
  } else if(!is.na(x) & is.na(y) & x<0){
    "Down regulated by IL4I1 and not differentially regulated in combination"
  } else if (is.na(x) & is.na(y)){
    "Check KD status"
  } else if(x > 0 & x > y){
    "Induced by IL4I1 and down regulated by KD in combination"
  } else if(x > 0 & x < y){
    "Induced by IL4I1 and even more by KD in combination"
  } else if (x < 0 & x < y){
    "Down regulated by IL4I1 and up regulated by KD in combination"
  } else if(x < 0 & x > y){
    "Down regulated by IL4I1 and even more by KD in combination"
  }
}))
SH2_comb_diff_apv <- unlist(map2(FCs_comp_apv$IL4I1_only, FCs_comp_apv$IL4I1_AHR_Dep_SH2,function(x, y){
  if(is.na(x) & !is.na(y) & y > 0){
    "Induced in combination _ check IL4I1 and KD contributions"
  } else if(is.na(x) & !is.na(y) & y < 0){
    "Down regulated in combination _ check IL4I1 and KD contributions"
  } else if (!is.na(x) & is.na(y) & x>0){
    "Induced by IL4I1  and not differentially regulated in combination"
  } else if(!is.na(x) & is.na(y) & x<0){
    "Down regulated by IL4I1 and not differentially regulated in combination"
  } else if (is.na(x) & is.na(y)){
    "Check KD status"
  } else if(x > 0 & x > y){
    "Induced by IL4I1 and down regulated by KD in combination"
  } else if(x > 0 & x < y){
    "Induced by IL4I1 and even more by KD in combination"
  } else if (x < 0 & x < y){
    "Down regulated by IL4I1 and up regulated by KD in combination"
  } else if(x < 0 & x > y){
    "Down regulated by IL4I1 and even more by KD in combination"
  }
}))

FCs_apv_df <- data.frame(FCs_comp_apv[,1:6], SH1_IL4I1_comp_group=SH1_comb_diff_apv,
                        SH2_IL4I1_comp_group=SH2_comb_diff_apv,
                        consensus_group = ifelse(SH1_comb_diff_apv == SH2_comb_diff_apv, "yes", "no"))

FCs_apv_df <- FCs_apv_df[-which(is.na(FCs_apv_df$IL4I1_only) & is.na(FCs_apv_df$IL4I1_AHR_Dep_SH1) & is.na(FCs_apv_df$SH1_no_IL4I1)),]
FCs_apv_df$IL4I1_Abs <- ifelse(!is.na(FCs_apv_df$IL4I1_only),ifelse(FCs_apv_df$IL4I1_only > 0, 2^FCs_apv_df$IL4I1_only, -1*(2^abs(FCs_apv_df$IL4I1_only))),NA)
FCs_apv_df$SH1_Abs <- ifelse(!is.na(FCs_apv_df$IL4I1_AHR_Dep_SH1),ifelse(FCs_apv_df$IL4I1_AHR_Dep_SH1 > 0, 2^FCs_apv_df$IL4I1_AHR_Dep_SH1, -1*(2^abs(FCs_apv_df$IL4I1_AHR_Dep_SH1))),NA)
FCs_apv_df$SH2_Abs <- ifelse(!is.na(FCs_apv_df$IL4I1_AHR_Dep_SH2),ifelse(FCs_apv_df$IL4I1_AHR_Dep_SH2 > 0, 2^FCs_apv_df$IL4I1_AHR_Dep_SH2, -1*(2^abs(FCs_apv_df$IL4I1_AHR_Dep_SH2))),NA)

FCs_apv_df$IL4I1_SH1_Abs_diff <- ifelse(!is.na(FCs_apv_df$IL4I1_Abs) & !is.na(FCs_apv_df$SH1_Abs),
                                       FCs_apv_df$IL4I1_Abs-FCs_apv_df$SH1_Abs, 
                                       ifelse(!is.na(FCs_apv_df$IL4I1_Abs) & is.na(FCs_apv_df$SH1_Abs),
                                              FCs_apv_df$IL4I1_Abs, 
                                              ifelse(is.na(FCs_apv_df$IL4I1_Abs) & !is.na(FCs_apv_df$SH1_Abs),
                                                     FCs_apv_df$SH1_Abs, NA)))
FCs_apv_df$IL4I1_SH2_Abs_diff <- ifelse(!is.na(FCs_apv_df$IL4I1_Abs) & !is.na(FCs_apv_df$SH2_Abs),
                                       FCs_apv_df$IL4I1_Abs-FCs_apv_df$SH2_Abs, 
                                       ifelse(!is.na(FCs_apv_df$IL4I1_Abs) & is.na(FCs_apv_df$SH2_Abs),
                                              FCs_apv_df$IL4I1_Abs, 
                                              ifelse(is.na(FCs_apv_df$IL4I1_Abs) & !is.na(FCs_apv_df$SH2_Abs),
                                                     FCs_apv_df$SH2_Abs, NA)))

FCs_apv_df$consensus_ABS <- ifelse(abs(FCs_apv_df$IL4I1_SH1_Abs_diff-FCs_apv_df$IL4I1_SH2_Abs_diff) <=0.5, "yes", "no")

FCs_apv_df$TFs <- names(TFSBDB_list_comb_genes)[match(FCs_apv_df$gene, names(TFSBDB_list_comb_genes))]
FCs_apv_df$AHR <- AHR_genes$Gene[match(FCs_apv_df$gene, AHR_genes$Gene)]

write.table(FCs_apv_df, "./TopTables/comparison_all_genes_grouping_apv.txt", row.names = F, quote = F, sep = "\t")

FCs_apv_df_AHR_only <- FCs_apv_df[!is.na(FCs_apv_df$AHR),]
write.table(FCs_apv_df_AHR_only, "./TopTables/comparison_all_genes_grouping_apv_AHR_only.txt", row.names = F, quote = F, sep = "\t")

##### Transcription Factor Analyses representations #####
## Not DONE ####
#### GSA ####
gsa_res <- gsa_FUN(top_t = tt_all_hgene_ord$IL4I1_only, rma_obj = rma_2sd, exprmnt = "IL4I1_shAHR", cntrst = "many",
                   cnt_mat = contrast_matrix, res_path = "./GSA/",wts = arr_wts,
                   msig.data.lists = msig.data.lists, msigs = 14, pltfrm = "Affy", d_m = design_eset)

for(i in 1:length(gsa_res)){
  walk(1:length(gsa_res[[i]]), safely(gsa_bar_plot_FUN), gsa_ls=gsa_res[[i]], res_path="/Figures/",
       coi=paste0("IL4I1_shAHR",names(gsa_res)[i]), wd=12, ht=8, pix=600)
}

##### Transcription Factor Analyses  enrichment scores #####
## Not DONE ####

#### AHR signature enrichment
roast_FUN <- function(tt_ids, cont_mat,gs,rma_obj, des_mat, arr_wts){
AHR_indexed_hi_lo_HPP <- ids2indices(gene.sets = gs, identifiers = tt_ids[,1])
prbsets_gsa_2sd_HPP <- rownames(tt_ids)
gsa_2sd_idx_HPP <- match(prbsets_gsa_2sd_HPP, rma_obj@featureData@data$probesetid)
gsa_2sd_eset_HPP <- rma_obj[gsa_2sd_idx_HPP]
AHR_hi_lo_roast_HPP <- roast(oligo::exprs(gsa_2sd_eset_HPP),weights=arr_wts,index = AHR_indexed_hi_lo_HPP, design = des_mat,
                             contrast = cont_mat, set.statistic = "floormean")
AHR_hi_lo_roast_HPP
}

all_roast <- map2(tts_2sd_u_AHR, as.data.frame(contrast_matrix), roast_FUN, gs=AHR_genes$Gene, rma_obj=rma_2sd, des_mat=design_eset,arr_wts=arr_wts) %>% do.call(rbind,.)
all_roast <- data.frame(Condition=rownames(all_roast), all_roast, stringsAsFactors = F)
write.table(all_roast, "./GSA/AHR_enrichment_AHR_KD_all.txt", row.names = F, sep = "\t")

pdf("./Figures/barcodeplot_m_t_stat_AHR_KD_all.pdf", width = 8, height = 12)
par(mfrow=c(3,2))
barcodeplot(tts_2sd_u_AHR$IL4I1_only$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$IL4I1_only$hGene)[[1]],main=all_roast$Condition[1])
barcodeplot(tts_2sd_u_AHR$IL4I1_AHR_Dep_SH1$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$IL4I1_AHR_Dep_SH1$hGene)[[1]],main=all_roast$Condition[2])
barcodeplot(tts_2sd_u_AHR$SH1_no_IL4I1$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$SH1_no_IL4I1$hGene)[[1]],main=all_roast$Condition[3])
barcodeplot(tts_2sd_u_AHR$IL4I1_AHR_Dep_SH2$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$IL4I1_AHR_Dep_SH2$hGene)[[1]],main=all_roast$Condition[4])
barcodeplot(tts_2sd_u_AHR$SH2_no_IL4I1$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR$SH2_no_IL4I1$hGene)[[1]],main=all_roast$Condition[5])
dev.off()



