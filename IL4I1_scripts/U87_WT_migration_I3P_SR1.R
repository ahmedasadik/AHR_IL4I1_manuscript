# This is script is for the analysis of the microarray of U87 wild type cells treated with I3P or SR1

library(AffyPipelineHG2ST)

data("TFSBDB_list","TFSBDB_list_comb","TFSBDB_list_comb_genes",
     "TFSBDB_list_genes","hgnc","msig.data.lists")

setwd("./IL4I1/Results_U87_migration")

dir.create("./RDS")
dir.create("./TopTables")
dir.create("./GSA")
dir.create("./Figures")
dir.create("./TFSDB")

#### generate raw_cel dataset ####
raw_path <- "/home/data/Raws/IL4I1/migration/"
list_cels <- list.files(raw_path, pattern = "CEL$", full.names = T)

raws_cel <- read.celfiles(list_cels)
saveRDS(raws_cel, "./RDS/raws_cel.rds")

#### Assigning sample names ####
s_names <- gsub("\\_\\(HuGene\\-2\\_0\\-st\\)\\.CEL","",sampleNames(raws_cel))
s_names <- gsub("I ","I_", s_names) %>% gsub(" ","",.) %>% gsub("-","_",.)
sampleNames(raws_cel) <- s_names

#### loading experimental condition covariates ####
exp_data_variables <- xlsx::read.xlsx("/home/data/Raws/IL4I1/Migration_20180117 U87 D, I3P 25 µM, I3P + SR1 1µM, SR1.xlsx",
                                      1, stringsAsFactors = F)
exp_data_variables <- exp_data_variables %>% .[-c(1,19),] %>% .[,-12]
colnames(exp_data_variables) <- exp_data_variables[1,]
exp_data_variables <- exp_data_variables[-1,]
exp_data_variables$Sample <-  gsub("I ","I_", exp_data_variables$Sample) %>% gsub(" ","",.) %>% gsub("-","_",.)
exp_data_variables <- exp_data_variables[match(s_names, exp_data_variables$Sample),]

#### Create phenodata & metadata ####
condition <- c(rep(c("ctrl", "I3P", "I3P_SR1", "SR1"),2), rep(c("ctrl", "I3P", "SR1", "I3P_SR1"), 2))
sample_id <- exp_data_variables$Sample
date_exp <- rep(c("d1", "d2", "d3", "d4"), each=4)
conc_rna <- exp_data_variables$`Conc (µg/µL)`
info <- data.frame(sample_id=sample_id, condition=condition, date_exp=date_exp, conc=conc_rna)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'date_exp', 'conc'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd

rma_2sd <- gen_rma_cel_nocross_median_FUN(raws_cel, sd_val = 2,snames = s_names,ct_off = 0.9)

### PCA plots to determine possible batch effects raws ####
pca_raw <- prcomp(t(oligo::exprs(rma(raws_cel))), scale. = T)
pdf("./Figures/pca_raws_Meth_Kyn.pdf", width = 6, height = 6)
factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA Migration groups", habillage=pd@data$condition, pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA Migration groups", habillage=pd@data$date_exp, pointsize = 6, mean.point = FALSE)

dev.off()

### PCA plots to determine possible batch effects ####
pca_rma <- prcomp(t(exprs(rma_2sd)), scale. = T)
pdf("./Figures/pca_normalized_Meth_Kyn.pdf", width = 6, height = 6)
factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA Migration groups", habillage=pd@data$condition, pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA Migration groups", habillage=pd@data$date_exp, pointsize = 6, mean.point = FALSE)
dev.off()

#### DGE ####
# Design matrix
design_eset <- model.matrix(~0+condition+date_exp, data=pd@data)
colnames(design_eset) <- gsub("condition","",colnames(design_eset))
design_eset <- cbind(design_eset, as.character(pd@data$conc)%>%as.numeric(.))
colnames(design_eset)[8] <- "conc"
arr_wts <- arrayWeights(oligo::exprs(rma_2sd), design = design_eset)
contrast_matrix <- makeContrasts(I3P-ctrl, SR1-ctrl, I3P_SR1-ctrl,  levels = design_eset)
## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset, weights = arr_wts)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_migration.rds")

## all tts
tt_2sd_all <- topTable(fit_2sd_eB, number = Inf,adjust.method = "BH")
tt_2sd_u_all <- tt_2sd_all[,c(29,25,40:45)]

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
walk2(tts_2sd_u, paste("./TopTables/tt_2sd_",names(tts_2sd_u),".txt",sep = ""), write.table, quote = F, sep = "\t")

## add AHR signature genes
AHR_genes <- read.delim("/home/analyses/Projects_analyses/AHR_IL4I1/AHR/Results/Signature/overlapping_AHR_signature.txt", sep = "\t")
tts_2sd_u_AHR <- map(tts_2sd_u, function(x,y){
  x$AHR_target <- match(x$hGene, y)
  x$AHR_target[!is.na(x$AHR_target)] <- "yes"
  x
}, y=AHR_genes$Gene)

walk2(tts_2sd_u_AHR, paste("./TopTables/tt_2sd_",names(tts_2sd_u),"_AHR.txt",sep = ""), write.table, quote = F, sep = "\t")

#### GSA ####
gsa_res <- gsa_FUN(top_t = tts_2sd_u_AHR$`I3P - ctrl`, rma_obj = rma_2sd, exprmnt = "migration", cntrst = "many",
                   cnt_mat = contrast_matrix, res_path = "./GSA/",wts = arr_wts,
                   msig.data.lists = msig.data.lists, msigs = 14, pltfrm = "Affy", d_m = design_eset)

for(i in 1:length(gsa_res)){
  walk(1:length(gsa_res[[i]]), safely(gsa_bar_plot_FUN), gsa_ls=gsa_res[[i]], res_path="/Figures/",
       coi=names(gsa_res)[i], wd=12, ht=8, pix=600)
}

roast_FUN <- function(tt_ids, cont_mat,gs,rma_obj, des_mat, arr_wts){
  AHR_indexed_hi_lo_HPP <- ids2indices(gene.sets = gs, identifiers = tt_ids[,1])
  prbsets_gsa_2sd_HPP <- rownames(tt_ids)
  gsa_2sd_idx_HPP <- match(prbsets_gsa_2sd_HPP, rma_obj@featureData@data$probesetid)
  gsa_2sd_eset_HPP <- rma_obj[gsa_2sd_idx_HPP]
  AHR_hi_lo_roast_HPP <- roast(oligo::exprs(gsa_2sd_eset_HPP),weights=arr_wts,index = AHR_indexed_hi_lo_HPP, design = des_mat,
                               contrast = cont_mat, set.statistic = "floormean")
  AHR_hi_lo_roast_HPP
}

all_roast <- map2(tts_2sd_u_AHR, as.data.frame(contrast_matrix), roast_FUN, gs=AHR_genes$Gene, rma_obj=rma_2sd, des_mat=design_eset, arr_wts=arr_wts) %>% do.call(rbind,.)
all_roast <- data.frame(Condition=rownames(all_roast), all_roast, stringsAsFactors = F)

write.table(all_roast, "./GSA/AHR_enrichment_migration_all.txt", row.names = F, sep = "\t")

pdf("./Figures/barcodeplot_m_t_stat_migration_all.pdf", width = 8, height = 12)
par(mfrow=c(3,1))
barcodeplot(tts_2sd_u_AHR[[1]]$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR[[1]]$hGene)[[1]],main=all_roast$Condition[1])
barcodeplot(tts_2sd_u_AHR[[2]]$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR[[2]]$hGene)[[1]],main=all_roast$Condition[2])
barcodeplot(tts_2sd_u_AHR[[3]]$t, ids2indices(gene.sets = AHR_genes$Gene, identifiers = tts_2sd_u_AHR[[3]]$hGene)[[1]],main=all_roast$Condition[3])
dev.off()



