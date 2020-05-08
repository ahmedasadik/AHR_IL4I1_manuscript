#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the comparison between metastatic and primary melanoma, showing higher IL4I1 and AHR activation in the metastatic melanoma.
#################################################

## Load libraries
library(purrr)
library(Hmisc)
library(edgeR)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)

## Load data
load("../Resources/TCGA-SKCMTranscriptome_ProfilingWed_Oct__4_161844_2017.RData")

## AHR signature
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt",stringsAsFactors = F)

## AHR activation comparison between the primary and metastatic melanomas
TCGA_SKCM_colData <- SummarizedExperiment::colData(data)
TCGA_SKCM_genes <- SummarizedExperiment::rowData(data)
TCGA_SKCM_counts <- SummarizedExperiment::assay(data)

samplesTM <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(TCGA_SKCM_counts),typesample = c("TM","TP"))
counts_TM <- TCGA_SKCM_counts[,match(samplesTM, colnames(TCGA_SKCM_counts))]
TCGA_SKCM_colData_TM <- TCGA_SKCM_colData[match(samplesTM, TCGA_SKCM_colData$barcode),]

## Create the DGE list
skcm_dge <- DGEList(counts_TM, samples = TCGA_SKCM_colData_TM, genes = TCGA_SKCM_genes)

###################################################
## IL4I1 level in primary vs metastatic patients ##
###################################################
## CPM the data
skcm_cpm <- cpm(skcm_dge, prior.count = 1, log = T)

## Create a data frame for plotting
met_pri_df <- data.frame(IL4I1=skcm_cpm[which(skcm_dge$genes$external_gene_name=="IL4I1"),],
                         pt=skcm_dge$samples$shortLetterCode, stringsAsFactors = F)
met_pri_df$pt <- factor(met_pri_df$pt, levels = c("TP","TM"))
pdf("../Results/Figures/IL4I1_met_vs_pr_melanoma.pdf", height = 8, width = 12)
ggboxplot(met_pri_df,"pt","IL4I1", xlab = FALSE, ylab = "IL4I1 (log2 CPM)", ylim=c(-3,12))+
  scale_y_continuous(name = "IL4I1 (log2 CPM)",limits = c(-3,12), breaks = c(0,4,8))+
  stat_compare_means(label.y = c(1.5,10))+ggbeeswarm::geom_quasirandom()
dev.off()

##################################################################
## Comparing AHR activation in primary and metastatic melanomas ##
##################################################################
## Create the design matrix
des_mat <- model.matrix(~0+as.factor(skcm_dge$samples$shortLetterCode))
colnames(des_mat) <- c("TM","TP")

## Filter set if not expressed in at least 20% of all samples and create voom object
A <- rowSums(cpm(skcm_dge) > 1) >=470*.2
dge <- skcm_dge[A,,keep.lib.size=FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, plot=FALSE)

## Contrast matrix
cont_mat <- makeContrasts(TM-TP, levels = des_mat)

## Differential regulation
vfit<-lmFit(v, des_mat)
vfit<-contrasts.fit(vfit,contrasts=cont_mat)
efit<-eBayes(vfit, robust = TRUE)
tt <- topTable(efit, number = Inf, adjust.method = "BH", sort.by = "M")

## Generate a barcodeplot
skcm_idx <- tt$external_gene_name %in% AHR_genes$Gene

pdf("../Results/Figures/barcodeplot_skcm_met_vs_pr.pdf", width = 8, height = 6, pointsize = 16)
barcodeplot(tt$t,skcm_idx)
dev.off()

## Run ROAST
v_idx_vector <- v$genes$external_gene_name %in% AHR_genes$Gene
skcm_comp_roast <- roast(v, index = v_idx_vector, design = des_mat,
                         contrast = cont_mat[,1], set.statistic = "floormean")
roast_cols <- c("NGenes",	"PropDown",	"PropUp",	"Direction",	"PValue",	"FDR",	"PValue.Mixed",	"FDR.Mixed")
AHR_roast_df <- as.matrix(rep(0,length(roast_cols)),nrow=1) %>% t() %>% as.data.frame()
colnames(AHR_roast_df) <- roast_cols
AHR_roast_df$NGenes <- skcm_comp_roast$ngenes.in.set
AHR_roast_df$PropDown <- skcm_comp_roast$p.value$Active.Prop[1] %>% round(.,4)
AHR_roast_df$PropUp <- skcm_comp_roast$p.value$Active.Prop[2] %>% round(.,4)
AHR_roast_df$Direction <- ifelse(AHR_roast_df$PropDown > AHR_roast_df$PropUp, "Down", "Up")
AHR_roast_df$PValue <- AHR_roast_df$FDR <- ifelse(AHR_roast_df$Direction=="Down",
                                                  round(skcm_comp_roast$p.value$P.Value[3], 4),
                                                  round(skcm_comp_roast$p.value$P.Value[3],4))
AHR_roast_df$PValue.Mixed <- AHR_roast_df$FDR.Mixed <- round(skcm_comp_roast$p.value$P.Value[4],4)
write.table(AHR_roast_df,"../Results/Tables/SKCM_met_vs_pr_roast.txt",sep = "\t")