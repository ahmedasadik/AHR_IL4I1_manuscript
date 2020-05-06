library(purrr)
library(ggplot2)
library(corrplot)
library(reshape)
library(gridExtra)
library(ggbeeswarm)
library(ggpubr)

## Read TCGA_TPMs
TCGA_TPMs <- readRDS("/home/data/Processed/TCGA/RNAseqV2/TCGA_TPMs.rds")

## estimate the mean or median of IDO1, IDO2, and TDO2 in the tumors and normal tissue
TCGA_enz_means <- map_df(TCGA_TPMs[-14], function(dge){
  IDO1_idx <- which(rownames(dge)=="IDO1")
  IDO2_idx <- which(rownames(dge)=="IDO2")
  TDO2_idx <- which(rownames(dge)=="TDO2")
  IL4I1_idx <- which(rownames(dge)=="IL4I1")
  tpm_m_sel <- dge[c(IDO1_idx, IDO2_idx, TDO2_idx, IL4I1_idx),] %>% t()
  colnames(tpm_m_sel) <- c("IDO1", "IDO2", "TDO2", "IL4I1")
  tpm_m_sel <- log2(tpm_m_sel)
  tpm_v_sel_med <- apply(tpm_m_sel, 2, median)
  names(tpm_v_sel_med) <- colnames(tpm_m_sel)
  tpm_v_sel_med
}) %>% as.data.frame() %>% t()

colnames(TCGA_enz_means) <- c("IDO1", "IDO2", "TDO2", "IL4I1")
TCGA_enz_means <- as.data.frame(TCGA_enz_means)
TCGA_enz_means[which(is.infinite(TCGA_enz_means$IDO2)),2] <- 0

rownames(TCGA_enz_means) <- gsub("TCGA_","", rownames(TCGA_enz_means))
pdf("/home/analyses/Projects_WF/General/BJC_review/corr_plot_Trp_enzymes_w_IL4I1_TCGA.pdf",width = 8, height = 12, pointsize = 24)
corrplot(corr = t(as.matrix(TCGA_enz_means)), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")
dev.off()

pdf("/home/analyses/Projects_WF/General/BJC_review/corr_plot_Trp_enzymes_wout_IL4I1_TCGA.pdf",width = 8, height = 12, pointsize = 24)
corrplot(corr = t(as.matrix(TCGA_enz_means[,-4])), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")
dev.off()

#################
## GTEX TPMs
gtex_tpms <- data.table::fread("/home/data/Processed/TCGA/TPMs/GTEX/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")
gtex_annot <- read.delim("/home/data/Raws/GTEX/annot", sep = "\t", stringsAsFactors = F)

GTEX_IDO1_idx <- which(gtex_tpms$Description=="IDO1")
GTEX_IDO2_idx <- which(gtex_tpms$Description=="IDO2")
GTEX_TDO2_idx <- which(gtex_tpms$Description=="TDO2")
GTEX_IL4I1_idx <- which(gtex_tpms$Description=="IL4I1")

gtex_tpms_enz <- gtex_tpms[c(GTEX_IDO1_idx, GTEX_IDO2_idx, GTEX_TDO2_idx, GTEX_IL4I1_idx),] %>% t()

colnames(gtex_tpms_enz) <- c("IDO1", "IDO2", "TDO2", "IL4I1")
gtex_tpms_enz <- gtex_tpms_enz[-c(1:2),]
gtex_tpms_enz_log2 <- apply(gtex_tpms_enz, 2, as.numeric) %>% log2()
rownames(gtex_tpms_enz_log2) <- rownames(gtex_tpms_enz)

new_gtex_annot <- gtex_annot[match(rownames(gtex_tpms_enz_log2), gtex_annot$SAMPID)%>%.[!is.na(.)],]

gtex_counts <- as.data.frame(gtex_tpms_enz_log2)[match(new_gtex_annot$SAMPID,rownames(gtex_tpms_enz_log2))%>%.[!is.na(.)],] %>% as.matrix()

gtex_counts_no_zs <- apply(gtex_counts[,1:4], 2, function(x){
  x[which(is.infinite(x))] <- 0
  x})

gtex_counts2 <- as.data.frame(gtex_counts_no_zs) %>% data.frame(., tissue= new_gtex_annot$SMTS)

gtex_tissues <- levels(as.factor(gtex_counts2$tissue))

gtex_to_plot <- map(gtex_tissues, function(gts, g_df){
  df <- as.matrix(g_df[which(g_df$tissue==gts),-5])
  df <- apply(df,2,as.numeric)
  tpm_v_sel_med <- apply(df, 2, median)
  names(tpm_v_sel_med) <- colnames(df)
  tpm_v_sel_med
}, g_df=gtex_counts2) %>% do.call("rbind",.)

rownames(gtex_to_plot) <- gtex_tissues

pdf("/home/analyses/Projects_WF/General/BJC_review/corr_plot_Trp_enzymes_w_IL4I1_GTEx.pdf",width = 8, height = 12, pointsize = 24)
corrplot(corr = t(as.matrix(gtex_to_plot)), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")
dev.off()

pdf("/home/analyses/Projects_WF/General/BJC_review/corr_plot_Trp_enzymes_wout_IL4I1_GTEx.pdf",width = 8, height = 12, pointsize = 24)
corrplot(corr = t(as.matrix(gtex_to_plot[,-4])), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")
dev.off()

## combined plot of enzyme expression in all tumors ####
pdf("/home/analyses/Projects_WF/General/BJC_review/corr_plot_Trp_enzymes_w_IL4I1_TCGA_GTEx_Comparison.pdf",width = 12, height = 6, pointsize = 24)
layout(c(1,2))
corrplot(corr = t(as.matrix(gtex_to_plot[,])), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")

corrplot(corr = t(as.matrix(TCGA_enz_means[,])), insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 2, cl.pos = "r",cl.cex = 0.5,
         is.corr = F, col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17), method = "circle")
dev.off()