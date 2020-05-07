# This is script is for the analysis of the microarray of U251 IL4I1 overexpressing cells

library(AffyPipelineHG2ST)
library(dplyr)

data("TFSBDB_list","TFSBDB_list_comb","TFSBDB_list_comb_genes",
     "TFSBDB_list_genes","hgnc","msig.data.lists")

setwd("./IL4I1/Results_U251_ovr_xprsn/")

dir.create("./RDS")
dir.create("./TopTables")
dir.create("./GSA")
dir.create("./Figures")
dir.create("./TFSDB")

#### generate raw_cel dataset ####
raw_path <- "/home/data/Raws/IL4I1/ctrl_ovrxprsn_U251/"
list_cels <- list.files(raw_path, pattern = "CEL$", full.names = T)
raws_cel <- read.celfiles(list_cels)
saveRDS(raws_cel, "./RDS/raws_cel_u251.rds")

#### Assigning sample names ####
s_names <- gsub("\\_\\(HuGene\\-2\\_0\\-st\\)\\.CEL","",sampleNames(raws_cel))
sampleNames(raws_cel) <- s_names

### loading experimental condition covariates ####
exp_data_variables <- read.csv("/home/data/Raws/IL4I1/20170822_U251 IL4I1 overexpressing cells_RNA microarray.csv",
                               sep = "\t", stringsAsFactors = F) %>% .[1:8,]

#### Create phenodata & metadata ####
condition <- rep(c("ctrl", "trt"), each =4)
sample_id <- c("c1", "c2", "c3","c4", "t1", "t2", "t3", "t4")
date_of_exp <- as.factor(exp_data_variables$Date.performed)
info <- data.frame(sample_id=sample_id, condition=condition, date=date_of_exp)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'date'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd
rma_2sd <- gen_rma_cel_nocross_median_FUN(raws_cel, sd_val = 2, snames = s_names, ct_off = 0.25)

### PCA plots to determine possible batch effects raws ####
pca_raw <- prcomp(t(exprs(rma(raws_cel))), scale. = T)
pdf("./Figures/pca_raws_u251.pdf", width = 6, height = 6)
factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$Date.performed))],
                 title = "PCA U251 date performed", habillage=exp_data_variables$Date.performed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$RNA))],
                 title = "PCA U251 RNA extraction", habillage=exp_data_variables$RNA,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA U251 samples", habillage=pd@data$condition,pointsize = 6, mean.point = FALSE)
dev.off()
### PCA plots to determine possible batch effects normalized ####
pca_rma <- prcomp(t(exprs(rma_2sd)), scale. = T)
pdf("./Figures/pca_normalized_u251.pdf", width = 6, height = 6)
factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$Date.performed))],
                 title = "PCA U251 date performed", habillage=exp_data_variables$Date.performed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$RNA))],
                 title = "PCA U251 RNA extraction", habillage=exp_data_variables$RNA,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA U251 samples", habillage=pd@data$condition,pointsize = 6, mean.point = FALSE)
dev.off()

#### Datasets without IL4I1 probes ####
# IL4I1 idcs
idx_IL4I1_2sd <- which(rma_2sd@featureData@data$geneid == "IL4I1")
# rma_cel without IL4I1
rma_2sd_IL4I1 <- rma_2sd[-idx_IL4I1_2sd]
saveRDS(rma_2sd_IL4I1, "./RDS/rma_2sd_IL4I1_u251.rds")

#### DGE ####
# Design matrix
design_eset <- model.matrix(~condition, data=pd@data)
colnames(design_eset)[2] <- c("ctrl_trt")

## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset)
fit_2sd_eB <- eBayes(fit_2sd, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_u251.rds")
tt_2sd <- topTable(fit_2sd_eB, coef=2, number = Inf, sort.by = "M", adjust.method = "BH")
tt_2sd_u <- tt_2sd[,c(29,25,40:45)]
write.table(tt_2sd_u, "./TopTables/IL4I1_tt_2sd_u251.txt", row.names = F, quote = F, sep = "\t")

## DGE 2sd -IL4I1
fit_2sd_IL4I1 <- lmFit(rma_2sd_IL4I1, design = design_eset)
fit_2sd_IL4I1_eB <- eBayes(fit_2sd_IL4I1, trend = T)
saveRDS(fit_2sd_IL4I1_eB, "./RDS/fit_2sd_IL4I1_eB_u251.rds")
tt_2sd_IL4I1 <- topTable(fit_2sd_IL4I1_eB, coef=2, number = Inf, sort.by = "M", adjust.method = "BH")
tt_2sd_IL4I1_u <- tt_2sd_IL4I1[,c(29,25,40:45)]
write.table(tt_2sd_IL4I1_u, "./TopTables/IL4I1_tt_2sd_before_annot_u251.txt", row.names = F, quote = F, sep = "\t")

## annotation update
tt_2sd_IL4I1_u <- annot_hgcn_FUN(tt_2sd_IL4I1_u, hgnc = hgnc)

## add AHR signature genes
AHR_genes <- read.delim("/home/analyses/Projects_analyses/AHR_IL4I1/AHR/Results/Signature/overlapping_AHR_signature.txt", sep = "\t")
tt_2sd_IL4I1_u$AHR_target <- match(tt_2sd_IL4I1_u$hGene, AHR_genes$Gene)
tt_2sd_IL4I1_u$AHR_target[!is.na(tt_2sd_IL4I1_u$AHR_target)] <- "yes"

write.table(tt_2sd_IL4I1_u, "./TopTables/IL4I1_tt_2sd_after_annot_u251.txt", row.names = F, quote = F, sep = "\t")

## Volcano_plot for AHR targets
vp_tt <- tt_2sd_IL4I1_u[!is.na(tt_2sd_IL4I1_u$AHR_target),]

pdf("./Figures/volcano_plot_top10_AHR_targets_U251.pdf", width = 12, height = 8, pointsize = 24)
ggplot(vp_tt, aes(x=logFC, y=(-log10(P.Value)),label=as.character(vp_tt$hGene)))+geom_point(size=5, color="red")+theme_bw()+
  theme(axis.text = element_text(size = 12),line = element_blank())+
  theme(legend.position="none")+ylab("-log10 p-value")+
  geom_text(aes(label=ifelse(vp_tt$logFC < -0.79,vp_tt$hGene,'')),hjust=-0.35, vjust=0.6, color="black")+
  geom_text(aes(label=ifelse(vp_tt$hGene %in% "HSPB2",vp_tt$hGene,'')),hjust=0.4, vjust=-0.75, color="black")+
  geom_text(aes(label=ifelse(vp_tt$hGene %in% c("DKK3", "LTBP1", "THBS1"),vp_tt$hGene,'')),hjust=-0.35, vjust=0.6, color="black")
dev.off()

##### Transcription Factor Analyses representations #####
no_cores <- detectCores() -2
cl <- makeCluster(no_cores)
clusterEvalQ(cl, {library(AffyPipelineHG2ST)})
clusterExport(cl,varlist = c("tt_2sd_IL4I1_u"))
TF_reps_comb <- parSapply(cl, TFSBDB_list_comb_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_IL4I1_u))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_IL4I1_u, rep_df=., fc_filt=0.58, res_path="/TFSDB/comb_")

TF_reps_comb$AHR_target <- match(TF_reps_comb$TFs, AHR_genes$Gene)
TF_reps_comb$AHR_target[!is.na(TF_reps_comb$AHR_target)] <- "yes"

write.table(TF_reps_comb, "./TFSDB/comb_comp_numbers_comparison_TFs_AHR_u251.txt", sep = "\t", row.names = F)

TF_reps_motifs <- parSapply(cl, TFSBDB_list_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_IL4I1_u))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_IL4I1_u, rep_df=., fc_filt=0.58, res_path="/TFSDB/motifs_")
write.table(TF_reps_motifs, "./TFSDB/motifs_comp_numbers_comparison_TFs_AHR_u251.txt", sep = "\t", row.names = F)
stopCluster(cl)

#### GSA ####
gsa_res <- gsa_FUN(tt_2sd_IL4I1_u, rma_2sd_IL4I1,d_m =design_eset,  exprmnt = "ctrl_ovrxprsn", cntrst = "simple", res_path = "./GSA/", msig.data.lists = msig.data.lists)
walk(1:length(gsa_res), gsa_bar_plot_FUN, gsa_ls=gsa_res, res_path="/Figures/",coi="IL4I1_ctrl_ovrxprsn_u251", wd=8, ht=6, pix=300)

##### Transcription Factor Analyses  enrichment scores #####
walk2(list(TFSBDB_list_comb, TFSBDB_list), c(FALSE, TRUE), function (a,b,c,d,e,f,g){
  TF_camera_gsa_FUN(TF_list = a,TF_motifs = b, top_T = c, rma_obj = d, cntrst = e, d_m = f, res_path = g)},
  c = tt_2sd_IL4I1_u, d = rma_2sd, e = 2, f = design_eset, g="./TFSDB/u251_")

## AHR_genes GSA using roast ####
AHR_indexed_hi_lo <- ids2indices(gene.sets = AHR_genes$Gene, identifiers = tt_2sd_IL4I1_u$hGene)
prbsets_gsa_2sd <- rownames(tt_2sd_IL4I1_u)
gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_2sd_IL4I1@featureData@data$probesetid)
gsa_2sd_eset <- rma_2sd_IL4I1[gsa_2sd_idx]

AHR_hi_lo_roast <- roast(exprs(gsa_2sd_eset),index = AHR_indexed_hi_lo, design = design_eset,
                         contrast = 2, set.statistic = "floormean")
write.table(AHR_hi_lo_roast, "./GSA/AHR_enrichment_U251.txt", sep = "\t")

pdf("./Figures/barcodeplot_m_t_stat_U251.pdf", width = 8, height = 6)
barcodeplot(tt_2sd_IL4I1_u$t, AHR_indexed_hi_lo$Set1)
dev.off()

## correlation matrix
rma_exprs <- exprs(rma_2sd)
exprs_idx <- which(rma_2sd@featureData@data$geneid == "IL4I1")
cor_matrix_sp <- Hmisc::rcorr(t(rma_exprs),rma_exprs[exprs_idx,], type = "spearman")
cor_matrix_pr <- Hmisc::rcorr(t(rma_exprs),rma_exprs[exprs_idx,])
saveRDS(cor_matrix_sp, "./cor_matrix_sp_ctrl_ovrxprsn_u251.rds")
saveRDS(cor_matrix_pr, "./cor_matrix_pr_ctrl_ovrxprsn_u251.rds")
## DONE
