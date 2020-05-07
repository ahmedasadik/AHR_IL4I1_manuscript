# This is script is for the analysis of the microarray of U87 wild type cells treated with I3P, HPP and PP

library(AffyPipelineHG2ST)
library(dplyr)

data("TFSBDB_list","TFSBDB_list_comb","TFSBDB_list_comb_genes",
     "TFSBDB_list_genes","hgnc","msig.data.lists")

setwd("./IL4I1/Results_U87_metabolites/")

dir.create("./RDS")
dir.create("./TopTables")
dir.create("./GSA")
dir.create("./Figures")
dir.create("./TFSDB")

#### generate raw_cel dataset ####
raw_path <- "/home/data/Raws/IL4I1/metabolites/"
list_cels <- list.files(raw_path, pattern = "CEL$", full.names = T)
raws_cel <- read.celfiles(list_cels)
saveRDS(raws_cel, "./RDS/raws_cel.rds")

#### Assigning sample names ####
s_names <- gsub("\\_\\(HuGene\\-2\\_0\\-st\\)\\.CEL","",sampleNames(raws_cel))
sampleNames(raws_cel) <- s_names

#### loading experimental condition covariates ####
exp_data_variables <- read.delim("/home/data/Raws/IL4I1/metabolites_metadata.txt",
                                 sep = ",", stringsAsFactors = F)
exp_data_variables <- exp_data_variables[order(exp_data_variables$Sample.name.Microarray),]

#### Create phenodata & metadata ####
condition <- rep(c("ctrl", "HPP", "I3P", "PP"), each =5)
sample_id <- c(paste0("ctrl_",1:5), paste0("HPP_",1:5), paste0("I3P_",1:5),paste0("PP_",1:5))
date_of_exp <- as.factor(exp_data_variables$Date.performed)
who_person <- exp_data_variables$Whow.did.it.
rna_person <- exp_data_variables$RNA
info <- data.frame(sample_id=sample_id, condition=condition, date=date_of_exp, who=who_person, rna=rna_person)
rownames(info) <- s_names
metadata <- data.frame(labelDescription=c('sample_id', 'condition', 'date', 'who', 'rna'), channel=factor('_ALL_'))
pd <- new('AnnotatedDataFrame', data=info, varMetadata=metadata)
phenoData(raws_cel) <- pd

rma_2sd <- gen_rma_cel_nocross_median_FUN(raws_cel, sd_val = 2,snames = s_names,ct_off = 0.25)

### PCA plots to determine possible batch effects raws ####
pca_raw <- prcomp(t(exprs(rma(raws_cel))), scale. = T)
pdf("./Figures/pca_raws_metabolites.pdf", width = 6, height = 6)
factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA metabolites conditions", habillage=pd@data$condition, pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$Date.performed))],
                 title = "PCA metabolites date performed", habillage=exp_data_variables$Date.performed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$Whow.did.it.))],
                 title = "PCA metabolites who performed it", habillage=exp_data_variables$Whow.did.it.,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_raw, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$RNA))],
                 title = "PCA metabolites who did RNA extraction", habillage=exp_data_variables$RNA,pointsize = 6, mean.point = FALSE)
dev.off()

### PCA plots to determine possible batch effects ####
pca_rma <- prcomp(t(exprs(rma_2sd)), scale. = T)
pdf("./Figures/pca_normalized_metabolites.pdf", width = 6, height = 6)
factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(pd@data$condition))],
                 title = "PCA metabolites conditions (nc_grouping)", habillage=pd@data$condition, pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$Date.performed))],
                 title = "PCA metabolites date performed (nc_grouping)", habillage=exp_data_variables$Date.performed,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$Whow.did.it.))],
                 title = "PCA metabolites who performed it (nc_grouping)", habillage=exp_data_variables$Whow.did.it.,pointsize = 6, mean.point = FALSE)

factoextra::fviz(pca_rma, element = "ind", geom = "point", color = brewer.pal(8, "Set1")[as.numeric(as.factor(exp_data_variables$RNA))],
                 title = "PCA metabolites who did RNA extraction (nc_grouping)", habillage=exp_data_variables$RNA,pointsize = 6, mean.point = FALSE)
dev.off()

#### DGE ####
# Design matrix
design_eset <- model.matrix(~0+condition+date, data=pd@data)
colnames(design_eset) <- gsub("condition","",colnames(design_eset))
contrast_matrix <- makeContrasts(HPP-ctrl, I3P-ctrl, PP-ctrl, levels = design_eset)
## DGE 2sd
fit_2sd <- lmFit(rma_2sd, design = design_eset)
fit_2sd_cont <- contrasts.fit(fit_2sd, contrasts = contrast_matrix)
fit_2sd_eB <- eBayes(fit_2sd_cont, trend = T)
saveRDS(fit_2sd_eB, "./RDS/fit_2sd_eB_metabolites_u87.rds")

## HPP tt
tt_2sd_HPP <- topTable(fit_2sd_eB, coef = 1,  sort.by = "M",  number = Inf,adjust.method = "BH")
tt_2sd_u_HPP <- tt_2sd_HPP[,c(29,25,40:45)]
write.table(tt_2sd_u_HPP, "./TopTables/metabolites_U87_tt_2sd_HPP.txt", row.names = F, quote = F, sep = "\t")
## annotation update
tt_2sd_u_HPP <- annot_hgcn_FUN(tt_2sd_u_HPP, hgnc = hgnc)
## add AHR signature genes
AHR_genes <- read.delim("/home/analyses/Projects_analyses/AHR_IL4I1/AHR/Results/Signature/overlapping_AHR_signature.txt", sep = "\t")
tt_2sd_u_HPP$AHR_target <- match(tt_2sd_u_HPP$hGene, AHR_genes$Gene)
tt_2sd_u_HPP$AHR_target[!is.na(tt_2sd_u_HPP$AHR_target)] <- "yes"
write.table(tt_2sd_u_HPP, "./TopTables/metabolites_U87_tt_2sd_HPP_after_annot.txt", row.names = F, quote = F, sep = "\t")

## I3P tt
tt_2sd_I3P <- topTable(fit_2sd_eB, coef = 2,  sort.by = "M",  number = Inf,adjust.method = "BH")
tt_2sd_u_I3P <- tt_2sd_I3P[,c(29,25,40:45)]
write.table(tt_2sd_u_I3P, "./TopTables/metabolites_U87_tt_2sd_I3P.txt", row.names = F, quote = F, sep = "\t")
## annotation update
tt_2sd_u_I3P <- annot_hgcn_FUN(tt_2sd_u_I3P, hgnc = hgnc)
## add AHR signature genes
tt_2sd_u_I3P$AHR_target <- match(tt_2sd_u_I3P$hGene, AHR_genes$Gene)
tt_2sd_u_I3P$AHR_target[!is.na(tt_2sd_u_I3P$AHR_target)] <- "yes"

write.table(tt_2sd_u_I3P, "./TopTables/metabolites_U87_tt_2sd_I3P_after_annot.txt", row.names = F, quote = F, sep = "\t")

## PP tt
tt_2sd_PP <- topTable(fit_2sd_eB, coef = 3,  sort.by = "M",  number = Inf,adjust.method = "BH")
tt_2sd_u_PP <- tt_2sd_PP[,c(29,25,40:45)]
write.table(tt_2sd_u_PP, "./TopTables/metabolites_U87_tt_2sd_PP.txt", row.names = F, quote = F, sep = "\t")
## annotation update
tt_2sd_u_PP <- annot_hgcn_FUN(tt_2sd_u_PP, hgnc = hgnc)
## add AHR signature genes
tt_2sd_u_PP$AHR_target <- match(tt_2sd_u_PP$hGene, AHR_genes$Gene)
tt_2sd_u_PP$AHR_target[!is.na(tt_2sd_u_PP$AHR_target)] <- "yes"

write.table(tt_2sd_u_PP, "./TopTables/metabolites_U87_tt_2sd_PP_after_annot.txt", row.names = F, quote = F, sep = "\t")

## Volcanoplots for comparison between the metabolites
vp_FUN <- function (tt_df, res_path, cndtn, lg_pos) {
  jpeg(filename = file.path(res_path, paste0("Volcano plot of ", cndtn, "_dge.jpg")), units = "in", width = 6, height = 6, res=300)
  # -log10(p.value)
  tt_df$logPV <- -log10(tt_df$P.Value)
  # assigne color code based on 0.01 pvalue and +/- 1.2 logFC
  tt_df$colorcode <- "grey80"
  tt_df$colorcode[which(tt_df$logFC <=-0.58 & tt_df$logPV >=1.3)] <- "grey55"
  tt_df$colorcode[which(tt_df$logFC >=0.58 & tt_df$logPV >=1.3)] <- "grey55"
  tt_df$colorcode[which(tt_df$AHR_target == "yes" & tt_df$logFC >=0.58 & tt_df$logPV >=1.3)] <- "red"
  tt_df$colorcode[which(tt_df$AHR_target == "yes" & tt_df$logFC <=-0.58 & tt_df$logPV >=1.3)] <- "red"
  # ideces of the up and down regulated AHR genes.
  up_idx <- which(tt_df$AHR_target == "yes" & tt_df$logFC > 0)
  dn_idx <- which(tt_df$AHR_target == "yes" & tt_df$logFC < 0)
  # plot volcano plot for the top 18 genes with logFC +/- 1.2
  plot(tt_df$logFC, tt_df$logPV, xlim=c(c(min(tt_df$logFC)-0.5), c(max(tt_df$logFC)+0.5)), ylim=c(0,c(max(tt_df$logPV)+0.5)), pch=19, col=tt_df$colorcode,
       xlab="log Fold Change", ylab="-log10(pvalue)", cex=0.9, main=paste0("Volcano plot of ", cndtn, "_dge"))
  abline(h = 1.3, v=c(-0.58, 0.58), col=c("grey"), lwd=2, lty=2)
  legend(lg_pos,bty = "n",legend = c("AHR targets"),
         col = c("red"),pch = 19, pt.cex = 0.9,cex = 0.8)
  dev.off()
}
vp_FUN(tt_df = tt_2sd_u_HPP, res_path = "./Figures/", cndtn = "HPP", lg_pos = "bottomright")
vp_FUN(tt_df = tt_2sd_u_I3P, res_path = "./Figures/", cndtn = "I3P", lg_pos = "bottomright")
vp_FUN(tt_df = tt_2sd_u_PP, res_path = "./Figures/", cndtn = "PP", lg_pos = "bottomright")

##### Transcription Factor Analyses representations #####
## HPP
no_cores <- detectCores() -2
cl <- makeCluster(no_cores)
clusterEvalQ(cl, {library(AffyPipelineHG2ST)})
clusterExport(cl,varlist = c("tt_2sd_u_HPP","tt_2sd_u_I3P","tt_2sd_u_PP"))
TF_reps_comb_HPP <- parSapply(cl, TFSBDB_list_comb_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_u_HPP))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_u_HPP, rep_df=., fc_filt=0.58, res_path="/TFSDB/metabolites_U87_HPP_comb_")

TF_reps_comb_HPP$AHR_target <- match(TF_reps_comb_HPP$TFs, AHR_genes$Gene)
TF_reps_comb_HPP$AHR_target[!is.na(TF_reps_comb_HPP$AHR_target)] <- "yes"

write.table(TF_reps_comb_HPP, "./TFSDB/metabolites_U87_HPP_comb_comp_numbers_comparison_TFs_AHR.txt", sep = "\t", row.names = F)

TF_reps_motifs_HPP <- parSapply(cl, TFSBDB_list_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_u_HPP))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_u_HPP, rep_df=., fc_filt=0.58, res_path="/TFSDB/metabolites_U87_HPP_motifs_")
write.table(TF_reps_motifs_HPP, "./TFSDB/metabolites_U87_HPP_motifs_comp_numbers_comparison_TFs_AHR.txt", sep = "\t", row.names = F)

## I3P
TF_reps_comb_I3P <- parSapply(cl, TFSBDB_list_comb_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_u_I3P))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_u_I3P, rep_df=., fc_filt=0.58, res_path="/TFSDB/metabolites_U87_I3P_comb_")

TF_reps_comb_I3P$AHR_target <- match(TF_reps_comb_I3P$TFs, AHR_genes$Gene)
TF_reps_comb_I3P$AHR_target[!is.na(TF_reps_comb_I3P$AHR_target)] <- "yes"

write.table(TF_reps_comb_I3P, "./TFSDB/metabolites_U87_I3P_comb_comp_numbers_comparison_TFs_AHR.txt", sep = "\t", row.names = F)

TF_reps_motifs_I3P <- parSapply(cl, TFSBDB_list_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_u_I3P))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_u_I3P, rep_df=., fc_filt=0.58, res_path="/TFSDB/metabolites_U87_I3P_motifs_")
write.table(TF_reps_motifs_I3P, "./TFSDB/metabolites_U87_I3P_motifs_comp_numbers_comparison_TFs_AHR.txt", sep = "\t", row.names = F)

## PP
TF_reps_comb_PP <- parSapply(cl, TFSBDB_list_comb_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_u_PP))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_u_PP, rep_df=., fc_filt=0.58, res_path="/TFSDB/metabolites_U87_PP_comb_")

TF_reps_comb_PP$AHR_target <- match(TF_reps_comb_PP$TFs, AHR_genes$Gene)
TF_reps_comb_PP$AHR_target[!is.na(TF_reps_comb_PP$AHR_target)] <- "yes"

write.table(TF_reps_comb_PP, "./TFSDB/metabolites_U87_PP_comb_comp_numbers_comparison_TFs_AHR.txt", sep = "\t", row.names = F)

TF_reps_motifs_PP <- parSapply(cl, TFSBDB_list_genes, function(x){
  unlist(invoke_map(list(TF_genes_tt_FUN),list(list(fc_filt=NA, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=NA, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=NA,rtrn_value = "dim"),
                                               list(fc_filt=0.58, pv=0.05, apv=0.2,rtrn_value = "dim")),
                    ptrns=x, t_table=tt_2sd_u_PP))}) %>% t() %>% TF_rep_final_df_FUN(tt_2sd_u_PP, rep_df=., fc_filt=0.58, res_path="/TFSDB/metabolites_U87_PP_motifs_")
write.table(TF_reps_motifs_PP, "./TFSDB/metabolites_U87_PP_motifs_comp_numbers_comparison_TFs_AHR.txt", sep = "\t", row.names = F)

stopCluster(cl)

#### GSA ####
gsa_res <- gsa_FUN(top_t = tt_2sd_u_HPP, rma_obj = rma_2sd, exprmnt = "metabolites_U87", cntrst = "many",
                   cnt_mat = contrast_matrix, res_path = "./GSA/",
                   msig.data.lists = msig.data.lists, msigs = 14, pltfrm = "Affy", d_m = design_eset)

for(i in 1:length(gsa_res)){
  walk(1:length(gsa_res[[i]]), safely(gsa_bar_plot_FUN), gsa_ls=gsa_res[[i]], res_path="/Figures/",
       coi=paste0("metabolites_U87",names(gsa_res)[i]), wd=12, ht=8, pix=600)
}

##### Transcription Factor Analyses  enrichment scores #####
for(i in 1:3){
  walk2(list(TFSBDB_list_comb, TFSBDB_list), c(FALSE, TRUE), safely(function (a,b,c,d,e,f,g){
    TF_camera_gsa_FUN(top_T = tt_2sd_u_HPP, TF_list = a, TF_motifs = b, pltfrm = c, rma_obj = d, cntrst = e, d_m = f, res_path = g)}),
    c="Affy", d = rma_2sd, e = contrast_matrix[,i], f = design_eset, g=paste0("./TFSDB/",colnames(contrast_matrix)[i],"_"))
}

#### AHR signature enrichment
## HPP
AHR_indexed_hi_lo_HPP <- ids2indices(gene.sets = AHR_genes$Gene, identifiers = tt_2sd_u_HPP$hGene)
prbsets_gsa_2sd_HPP <- rownames(tt_2sd_u_HPP)
gsa_2sd_idx_HPP <- match(prbsets_gsa_2sd_HPP, rma_2sd@featureData@data$probesetid)
gsa_2sd_eset_HPP <- rma_2sd[gsa_2sd_idx_HPP]

AHR_hi_lo_roast_HPP <- roast(exprs(gsa_2sd_eset_HPP),index = AHR_indexed_hi_lo_HPP, design = design_eset,
                         contrast = contrast_matrix[,1], set.statistic = "floormean")
write.table(AHR_hi_lo_roast_HPP, "./GSA/AHR_enrichment_U87_HPP.txt", sep = "\t")

pdf("./Figures/barcodeplot_m_t_stat_U87_HPP.pdf", width = 8, height = 6)
barcodeplot(tt_2sd_u_HPP$t, AHR_indexed_hi_lo_HPP$Set1)
dev.off()

## I3P
AHR_indexed_hi_lo_I3P <- ids2indices(gene.sets = AHR_genes$Gene, identifiers = tt_2sd_u_I3P$hGene)
prbsets_gsa_2sd_I3P <- rownames(tt_2sd_u_I3P)
gsa_2sd_idx_I3P <- match(prbsets_gsa_2sd_I3P, rma_2sd@featureData@data$probesetid)
gsa_2sd_eset_I3P <- rma_2sd[gsa_2sd_idx_I3P]

AHR_hi_lo_roast_I3P <- roast(exprs(gsa_2sd_eset_I3P),index = AHR_indexed_hi_lo_I3P, design = design_eset,
                             contrast = contrast_matrix[,2], set.statistic = "floormean")
write.table(AHR_hi_lo_roast_I3P, "./GSA/AHR_enrichment_U87_I3P.txt", sep = "\t")

pdf("./Figures/barcodeplot_m_t_stat_U87_I3P.pdf", width = 8, height = 6)
barcodeplot(tt_2sd_u_I3P$t, AHR_indexed_hi_lo_I3P$Set1)
dev.off()

## PP
AHR_indexed_hi_lo_PP <- ids2indices(gene.sets = AHR_genes$Gene, identifiers = tt_2sd_u_PP$hGene)
prbsets_gsa_2sd_PP <- rownames(tt_2sd_u_PP)
gsa_2sd_idx_PP <- match(prbsets_gsa_2sd_PP, rma_2sd@featureData@data$probesetid)
gsa_2sd_eset_PP <- rma_2sd[gsa_2sd_idx_PP]

AHR_hi_lo_roast_PP <- roast(exprs(gsa_2sd_eset_PP),index = AHR_indexed_hi_lo_PP, design = design_eset,
                             contrast = contrast_matrix[,3], set.statistic = "floormean")
write.table(AHR_hi_lo_roast_PP, "./GSA/AHR_enrichment_U87_PP.txt", sep = "\t")

pdf("./Figures/barcodeplot_m_t_stat_U87_PP.pdf", width = 8, height = 6)
barcodeplot(tt_2sd_u_PP$t, AHR_indexed_hi_lo_PP$Set1)
dev.off()

## combined roast results
roast_all <- rbind(AHR_hi_lo_roast_HPP, AHR_hi_lo_roast_I3P, AHR_hi_lo_roast_PP)
rownames(roast_all) <- c("HPP", "I3P", "PP")
write.table(roast_all, "./GSA/AHR_enrichment_U87_all_mets.txt", sep = "\t")


