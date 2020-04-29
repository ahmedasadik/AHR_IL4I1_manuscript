# This script describes the generation of toptables of all comparisons for k=2:5 across all tumors after CC

## Load libraries
library(purrr)
library(edgeR)

# Concensus clustering results
ccp_pos_dist <- readRDS("./RDS/ccp_pos_dist.rds")

# TCGA_voom
TCGA_voom <- readRDS("./RDS/TCGA_DGE_voom_annots.rds")

# The selected number of clusters of each TCGA tumor (this was analyzed manually as previously described)
TCGA_clusters <- paste("K",c(2,4,4,4,4,4,2,4,4,4,2,3,3,4,3,4,4,2,4,2,2,2,4,4,2,2,3,3,3,3,2,3), sep = "")

## Differential gene expression function
RNA_DGE_WGCNA_voomed_data_FUN <- function(des_mat, dge){
  if(length(colnames(des_mat))==2){
    cont_mat <- makeContrasts(cls_1-cls_2, levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit, robust = TRUE)
    tt <- topTable(efit, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt, eBFit=efit)
  } else if(length(colnames(des_mat))==3){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_2-cls_3, levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit)
    tt1 <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "M")
    tt2 <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")
    tt3 <- topTable(efit, coef = 3, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt1, cnt2=tt2, cnt3=tt3, eBFit=efit)
  } else if(length(colnames(des_mat))==4){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_1-cls_4, cls_2-cls_3, cls_2-cls_4, cls_3-cls_4,levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit)
    tt1 <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "M")
    tt2 <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")
    tt3 <- topTable(efit, coef = 3, number = Inf, adjust.method = "BH", sort.by = "M")
    tt4 <- topTable(efit, coef = 4, number = Inf, adjust.method = "BH", sort.by = "M")
    tt5 <- topTable(efit, coef = 5, number = Inf, adjust.method = "BH", sort.by = "M")
    tt6 <- topTable(efit, coef = 6, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt1, cnt2=tt2, cnt3=tt3, cnt4=tt4, cnt5=tt5, cnt6=tt6, eBFit=efit)
  } else if(length(colnames(des_mat))==5){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_1-cls_4,cls_1-cls_5,
                              cls_2-cls_3, cls_2-cls_4,cls_2-cls_5, cls_3-cls_4,
                              cls_3-cls_5,cls_4-cls_5,levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit)
    tt1 <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "M")
    tt2 <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")
    tt3 <- topTable(efit, coef = 3, number = Inf, adjust.method = "BH", sort.by = "M")
    tt4 <- topTable(efit, coef = 4, number = Inf, adjust.method = "BH", sort.by = "M")
    tt5 <- topTable(efit, coef = 5, number = Inf, adjust.method = "BH", sort.by = "M")
    tt6 <- topTable(efit, coef = 6, number = Inf, adjust.method = "BH", sort.by = "M")
    tt7 <- topTable(efit, coef = 7, number = Inf, adjust.method = "BH", sort.by = "M")
    tt8 <- topTable(efit, coef = 8, number = Inf, adjust.method = "BH", sort.by = "M")
    tt9 <- topTable(efit, coef = 9, number = Inf, adjust.method = "BH", sort.by = "M")
    tt10 <- topTable(efit, coef =10, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt1, cnt2=tt2, cnt3=tt3, cnt4=tt4, cnt5=tt5, cnt6=tt6,
                cnt7=tt7, cnt8=tt8, cnt9=tt9, cnt10=tt10,
                eBFit=efit)
  }
}

## Differential gene expression using using the voomed data sets and the CC results.
## The design and contrast matrices are automatically inferred
TCGA_DGE_ccp <- purrr::map2(TCGA_voom, ccp_pos_dist, function(a1, a2){
  ccp_input <- a2$ccp
  ccp_idcs <- cbind(2,3,4,5)
  ccp_idcs_mat <- apply(ccp_idcs,2, function(x){
    ccp_input[[x]]$consensusClass
  })
  des_mats <- apply(ccp_idcs_mat,2,function(y){
    des<-model.matrix(~0+as.factor(y))
    colnames(des) <- paste("cls",levels(as.factor(y)),sep="_")
    des})
  names(des_mats) <- c("K2", "K3", "K4", "K5")
  tts <- lapply(des_mats, RNA_DGE_WGCNA_voomed_data_FUN, dge=a1)
})

saveRDS(TCGA_DGE_ccp, "./RDS/TCGA_DGE_ccp.rds")

## selected TCGA_clusters comparisons
TCGA_DGE_ccp_clusters <- map2(TCGA_DGE_ccp, TCGA_clusters, function(a1,a2){
  tts <- a1[[a2]] %>% .[-length(.)]
})

saveRDS(TCGA_DGE_ccp_clusters, "./RDS/TCGA_DGE_ccp_clusters.rds")

walk2(TCGA_DGE_ccp_clusters, names(TCGA_DGE_ccp_clusters), function(a1,a2){
  walk2(a1, paste("./Tables/",a2,"_tt_", names(a1),".csv",sep = ""), write.csv)
})
