#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript/AHR_signature/
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This files describes how the AHR signature was generated
################################################

## Load libraries
library(limma)
library(purrr)
library(lumi)
library(GEOquery)
library(Biobase)
library(ggplot2)


## AHR immune mouse signature
AHR_mouse_immune <- read.delim("./AHR_mouse_immune_GRCm38.txt", sep = "\t", stringsAsFactors = F)

## preprocessing ####
sample_sheet <- read.delim("./sample_sheet_hannab_17.03.2014.csv",skip = 11, sep = ",") %>% .[,-c(2,3,5)] %>% .[1:10,]

all_raws <- read.csv("./all.raw.csv", sep = "\t", stringsAsFactors = F)
ILMN_raws <- all_raws[,c(65, 3:42, 52:58, 60:64,70, 74,78)]
cll_mc_raws <- ILMN_raws[,c(1, 52, 2:41, 44,46,54,55,56)]

colnames(cll_mc_raws) <- toupper(colnames(cll_mc_raws))
colnames(cll_mc_raws)[3:42] <- unlist(map(1:10,function(x){paste(x, c("AVG_Signal", "Avg_NBEADS", "BEAD_STDERR", "Detection.Pval"), sep = ".")}))
rownames(cll_mc_raws) <- cll_mc_raws$PROBE_ID

write.table(cll_mc_raws, "../Results/GSE66858/spleen_monocytes_raws.txt", sep = "\t", row.names = F)

gse_path <- "../Results/GSE66858/spleen_monocytes_raws.txt"
lumi_obj <- lumiR(gse_path)
saveRDS(lumi_obj, "../Results/GSE66858/lumi_obj_spleen_mice_monocytes_total.rds")

lumi_obj_T <- lumiT(lumi_obj, "vst")
lumi_obj_N <- lumiN(lumi_obj_T,method = "rsn")
lumi.N.Q <- lumiQ(lumi_obj_N)
saveRDS(lumi.N.Q, "../Results/GSE66858/lumi_N.Q_spleen_mice_monocytes_total.rds")

## filtering unexpressed probes
expressed <- detectionCall(lumi.N.Q)
edata <- exprs(lumi.N.Q)
edata <- edata[expressed > 0,]
saveRDS(edata, "../Results/GSE66858/edata_spleen_mice_monocytes_total.rds")

## Annotating probes
## matching annotated probes to lumi object
map_idx <- match(rownames(edata), lumi.N.Q@featureData@data$PROBE_ID)
rownames(edata) <- lumi.N.Q@featureData@data$SYMBOL[map_idx]

## expression matrix
edata_av <- avereps(edata)

saveRDS(edata_av, "../Results/GSE66858/edata_av_spleen_mice_monocytes_total.rds")

## Differential expression limma ####
ct <- factor(sample_sheet$Sample_Group)
ct <- relevel(ct, ref="WT")
design_mat <- model.matrix(~ct)
colnames(design_mat) <- gsub("ct", "", colnames(design_mat))
fit <- lmFit(edata_av,design_mat)
fit2 <- eBayes(fit, trend=TRUE)
tt <- topTable(fit2, coef=2, number = Inf,adjust.method = "BH", sort.by = "M")
tt <- data.frame(hGene=rownames(tt), tt, stringsAsFactors = F)
tt <- tt [order(tt$logFC, decreasing = T),]
## Volcano plot ##

## Volcano_plot for AHR targets
vp_tt <- tt
vp_tt$cc <- "grey"
vp_tt$cc[which(vp_tt$logFC > 2.5)] <- "red"

pdf("../Results/GSE66858/GSE66858_volcanoplot.pdf",width = 8, height = 6, pointsize = 24)
ggplot(vp_tt, aes(x=logFC, y=(-log10(P.Value)),label=as.character(vp_tt$hGene)))+geom_point(size=5, color=vp_tt$cc)+theme_bw()+
  theme(axis.text = element_text(size = 12),line = element_blank())+
  theme(legend.position="none")+ylab("-log10 p-value")+
  geom_vline(xintercept = 0.58, linetype="dashed", color = "grey")+
  geom_vline(xintercept = -0.58, linetype="dashed", color = "grey")+
  geom_hline(yintercept=2, linetype="dashed", color = "grey")+
  geom_text(aes(label=ifelse((vp_tt$hGene %in% vp_tt$hGene[c(1:3)]),vp_tt$hGene,'')),hjust=1.1, vjust=-0.5, color="black")+
  geom_text(aes(label=ifelse((vp_tt$hGene %in% vp_tt$hGene[4]),vp_tt$hGene,'')),hjust=-0.3, vjust=-0.5, color="black")+
  geom_text(aes(label=ifelse((vp_tt$hGene %in% vp_tt$hGene[c(5:7)]),vp_tt$hGene,'')),hjust=1.1, vjust=-0.5, color="black")
dev.off()

## Barcodeplots
AHR_indexed_immune <- tt$hGene %in% AHR_mouse_immune$hGene
pdf("../Results/GSE66858/GSE66858_barcodeplot_mouse_CLL.pdf", width = 8, height = 6)
barcodeplot(tt$t, AHR_indexed_immune, main="GSE66858_CLL_Monocytes")
dev.off()

## ROAST
prbsets_immune <- rownames(tt)
gsa_immune <- match(prbsets_immune, rownames(edata_av))
gsa_bc_immune <- edata_av[gsa_immune,]
AHR_bc_roast_immune <- roast(gsa_bc_immune, index = AHR_indexed_immune, design = design_mat, contrast = 2, set.statistic = "floormean")
roast_cols <- c("NGenes",	"PropDown",	"PropUp",	"Direction",	"PValue",	"FDR",	"PValue.Mixed",	"FDR.Mixed")
AHR_roast_df <- as.matrix(rep(0,length(roast_cols)),nrow=1) %>% t() %>% as.data.frame()
colnames(AHR_roast_df) <- roast_cols
AHR_roast_df$NGenes <- AHR_bc_roast_immune$ngenes.in.set
AHR_roast_df$PropDown <- AHR_bc_roast_immune$p.value$Active.Prop[1] %>% round(.,4)
AHR_roast_df$PropUp <- AHR_bc_roast_immune$p.value$Active.Prop[2] %>% round(.,4)
AHR_roast_df$Direction <- ifelse(AHR_roast_df$PropDown > AHR_roast_df$PropUp, "Down", "Up")
AHR_roast_df$PValue <- AHR_roast_df$FDR <- ifelse(AHR_roast_df$Direction=="Down",
                                                  round(AHR_bc_roast_immune$p.value$P.Value[3], 4),
                                                  round(AHR_bc_roast_immune$p.value$P.Value[3],4))
AHR_roast_df$PValue.Mixed <- AHR_roast_df$FDR.Mixed <- round(AHR_bc_roast_immune$p.value$P.Value[4],4)

write.table(AHR_roast_df, "../Results/GSE66858/GSE66858_roast_res.txt")
