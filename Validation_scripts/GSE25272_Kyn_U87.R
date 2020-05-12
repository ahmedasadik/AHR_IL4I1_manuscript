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

library(limma)
library(purrr)
library(lumi)
library(GEOquery)
library(Biobase)
library(illuminaHumanv3.db)

## Load the AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## Load the raw data
all_raws <- read.delim("../Validation_datasets/GSE25272/GSE25272_non-normalized.txt", sep = "\t", stringsAsFactors = F)

## preprocessing ####
Sentrix_Position <- LETTERS[1:8]
Sentrix_ID <- rep(440382705,8)
Sample_Name <- c("ut_8_1","ut_8_2","kyn_8_1","kyn_8_2","ut_24_1","ut_24_2","kyn_24_1","kyn_24_2")
Sample_group <- rep(c("c8","k8","c24","k24"),each=2)
sample_sheet <- data.frame(Sample_Name, Sample_group, Sentrix_ID, Sentrix_Position, stringsAsFactors = F)

ILMN_raws <- all_raws[-c(which(nchar(all_raws$Transcript)==0)),]
ILMN_raws_ordered <- all_raws[,c(1, 40, 3:34, 38,40,54,58,62)]

colnames(ILMN_raws_ordered) <- toupper(colnames(ILMN_raws_ordered))
colnames(ILMN_raws_ordered)[3:34] <- unlist(map(1:8,function(x){paste(x, c("AVG_Signal", "BEAD_STDV",  "Detection.Pval","RELATIVE_SD"), sep = ".")}))
colnames(ILMN_raws_ordered)[1] <- "PROBE_ID"
rownames(ILMN_raws_ordered) <- ILMN_raws_ordered$PROBE_ID

write.table(ILMN_raws_ordered, "../Results/Validations/GSE25272_U87/GSE25272_ILMN_raws.txt", sep = "\t", row.names = F)

gse_path <- "../Results/Validations/GSE25272_U87/GSE25272_ILMN_raws.txt"
lumi_obj <- lumiR(gse_path)
lumi_obj_T <- lumiT(lumi_obj, "vst")
lumi_obj_N <- lumiN(lumi_obj_T,method = "rsn")
lumi.N.Q <- lumiQ(lumi_obj_N)
saveRDS(lumi.N.Q, "../Results/Validations/GSE25272_U87/GSE25272_lumi_obj.rds")

## filtering unexpressed probes
expressed <- detectionCall(lumi.N.Q)
edata <- exprs(lumi.N.Q)
edata <- edata[expressed > 0,]

## Annotating probes
## matching annotated probes to lumi object
probes_w_transcripts <- all_raws[-c(which(nchar(all_raws$Transcript)==0)),c("ProbeID", "ILMN_Gene", "Transcript")] %>% data.frame(.,stringsAsFactors=F)
map_idx <- rownames(edata) %in% probes_w_transcripts$ProbeID

edata <- edata[map_idx,]

gene_idx <- probes_w_transcripts$ProbeID %in% rownames(edata)
rownames(edata) <- probes_w_transcripts$ILMN_Gene[gene_idx]

saveRDS(edata, "../Results/Validations/GSE25272_U87/GSE25272_edata_obj.rds")

## expression matrix
edata_av <- avereps(edata)
colnames(edata_av) <- sample_sheet$Sample_Name
saveRDS(edata_av, "../Results/Validations/GSE25272_U87/GSE25272_edata_av_obj.rds")

## Run DGE
ct <- factor(sample_sheet$Sample_group)
design_mat <- model.matrix(~0+ct)
colnames(design_mat) <- gsub("ct", "", colnames(design_mat))
cont_mat <- makeContrasts(k8-c8, k24-c24, levels = design_mat)
fit <- lmFit(edata_av,design_mat)
fit <- contrasts.fit(fit, contrasts = cont_mat)
fit2 <- eBayes(fit, trend=TRUE, robust = TRUE)
tt_all <- topTable(fit2, number = Inf,adjust.method = "BH")

tt_all_8 <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", sort.by = "M")
tt_all_24 <- topTable(fit2, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")

## Barcodeplot
index_all_8 <- rownames(tt_all_8) %in% AHR_genes$Gene
pdf("../Results/Validations/GSE25272_U87/barcodeplot_U87.pdf",width = 12, height = 8)
barcodeplot(tt_all_8$t, index_all_8, main="GSE25272_Kyn_U87_8H")
dev.off()

ilmn_roast_FUN <- function(v_eset, glist, des_mat, cont_mat, ...){
  index.vector <- rownames(v_eset) %in% glist
  AHR_roast <- roast(v_eset, index = index.vector, design = des_mat,
                     contrast = cont_mat, set.statistic = "floormean")
  roast_cols <- c("NGenes",	"PropDown",	"PropUp",	"Direction",	"PValue",	"FDR",	"PValue.Mixed",	"FDR.Mixed")
  AHR_roast_df <- as.matrix(rep(0,length(roast_cols)),nrow=1) %>% t() %>% as.data.frame()
  colnames(AHR_roast_df) <- roast_cols
  AHR_roast_df$NGenes <- AHR_roast$ngenes.in.set
  AHR_roast_df$PropDown <- AHR_roast$p.value$Active.Prop[1] %>% round(.,4)
  AHR_roast_df$PropUp <- AHR_roast$p.value$Active.Prop[2] %>% round(.,4)
  AHR_roast_df$Direction <- ifelse(AHR_roast_df$PropDown > AHR_roast_df$PropUp, "Down", "Up")
  AHR_roast_df$PValue <- AHR_roast_df$FDR <- ifelse(AHR_roast_df$Direction=="Down",
                                                    round(AHR_roast$p.value$P.Value[3], 4),
                                                    round(AHR_roast$p.value$P.Value[3],4))
  AHR_roast_df$PValue.Mixed <- AHR_roast_df$FDR.Mixed <- round(AHR_roast$p.value$P.Value[4],4)
  AHR_roast_df 
}

cont1_roast <- ilmn_roast_FUN(edata_av, AHR_genes$Gene, design_mat, cont_mat[,1], trend=TRUE)
cont2_roast <- ilmn_roast_FUN(edata_av, AHR_genes$Gene, design_mat, cont_mat[,2], trend=TRUE)

AHR_roast_df <- data.frame(Condition=c("8H","24H"), rbind(cont1_roast,cont2_roast),stringsAsFactors = F)

write.table(AHR_roast_df, "../Results/Validations/GSE25272_U87/GSE25272_roast_res.txt", sep = "\t")
