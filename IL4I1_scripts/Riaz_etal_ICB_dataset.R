#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes how we came across IL4I1 as an enzyme producing ligands activating AHR.
#################################################

library(edgeR)
library(purrr)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)

## Source the functions and parameters files
source("../functions_and_parameters.R")

## The files for this dataset were downloaded from the publication github repo on october 2018 
## (https://github.com/riazn/bms038_analysis/tree/master/)

## read raw counts data
riaz_counts <- read.delim("../Resources/CountData.BMS038.txt")
riaz_genes <- data.frame(genes=riaz_counts$HUGO, eids=rownames(riaz_counts), stringsAsFactors = F)
riaz_counts_numeric <- apply(riaz_counts[,-119],2,as.numeric)
AHR_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

## read pheno data
riaz_samples <- read.csv("../Resources/SampleTableCorrected.9.19.16.csv", stringsAsFactors = F)
riaz_samples <- riaz_samples[,-1]
riaz_clinical <- read.csv("../Resources/bms038_clinical_data.csv", stringsAsFactors = F)

## RIAZ DGElist
riaz_counts_numeric <- riaz_counts_numeric[,colnames(riaz_counts_numeric) %in% riaz_samples$Sample]
riaz_samples <- riaz_samples[riaz_samples$Sample%in% colnames(riaz_counts_numeric),]

riaz_samples$subtype <- ""
for(i in 1:length(riaz_samples$subtype)){
  riaz_samples$subtype[i] <- riaz_clinical$SubtypeEZ[match(riaz_samples$PatientID[i], riaz_clinical$PatientID)]
}
dge <- DGEList(counts = riaz_counts_numeric, genes = riaz_genes, samples = riaz_samples)

## We are using the filtered and normalized RNA_seq data to reduce the effect of noise
keep <- filterByExpr(dge)
dge_voom <- dge[keep,,keep.lib.size=FALSE]
dge_voom <- calcNormFactors(dge_voom)
dge_voom <- voom(dge_voom, plot=FALSE)
rownames(dge_voom$E) <- dge_voom$genes$genes

dge_voom$targets$PreOn <- factor(dge_voom$targets$PreOn, levels = c("Pre", "On"))

## Split the patients into different groups based on the treatment plan, subtype and response
dge_voom$targets$group_ttt_subtype <- paste(dge_voom$targets$Cohort, dge_voom$targets$subtype, sep = "_")
dge_voom$targets$group_ttt_subtype_response <- paste(dge_voom$targets$Cohort, dge_voom$targets$subtype, dge_voom$targets$BOR, sep = "_")

## Perform DGE simple before and after treatment  All patients ####
des_mat <- model.matrix(~PreOn, dge_voom$targets)

#des_mat <- model.matrix(~PatientID+PreOn, dge_voom$targets)
corfit <- duplicateCorrelation(dge_voom, design = des_mat,block=dge_voom$targets$PatientID)

## toptable
vfit<-lmFit(dge_voom, des_mat, block=dge_voom$targets$PatientID, correlation = corfit$consensus)
efit <- eBayes(vfit)
tt <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")
tt <- tt[order(tt$logFC, decreasing = T),]

write.csv(tt,"../Results/Tables/RIAZ_ICB_dataset_tt.csv", row.names = F)

## The threshold of significance is 0.01
tt_plot <- tt[tt$genes %in% c("IDO1","IL4I1", "TDO2"),]
tt_plot$FC <- 2^tt_plot$logFC %>% round(.,1)
tt_plot$FC <- ifelse(tt_plot$P.Value<0.01, tt_plot$FC, 0.001)
tt_plot$fill_colors <- ifelse(tt_plot$FC > 0.1, "1","2")

tt_plot$genes <- factor(tt_plot$genes, levels = c("IDO1","IL4I1", "TDO2"))

tt_p <- ggbarplot(tt_plot, x = "genes", y="FC", fill = "fill_colors", palette = "npg", xlab = FALSE,
                  ylab = "Fold Change\n", legend.position="none")+theme(text = element_text(size = 12))

label.df <- data.frame(genes = c("IDO1", "IL4I1", "TDO2"), FC = c(2.2, 1.7, 0.3))
tt_p2 <- tt_p + geom_text(data = label.df, label = c("****","***","ns"), size=6)+theme(legend.position = "none")

pdf("../Results/Figures/Riaz_ICB_dataset_barplot_IL4I1_IDO1_TDO2.pdf",width =3, height = 4,pointsize = 12)
tt_p2
dev.off()

tt_idx <- tt$genes %in% AHR_genes$Gene

pdf("../Results/Figures/Riaz_ICB_dataset_barcode_plot.pdf", width = 12, height = 8, pointsize = 16)
barcodeplot(tt$t, tt_idx, font=5, xlab = "")
dev.off()

index.vector <- rownames(dge_voom$E) %in% AHR_genes$Gene
AHR_roast <- roast(dge_voom, index = index.vector, design = des_mat, block=dge_voom$targets$PatientID,
                   correlation = corfit$consensus,
                   contrast = 2)
roast_cols <- c("NGenes",	"PropDown",	"PropUp",	"Direction",	"PValue",	"FDR",	"PValue.Mixed",	"FDR.Mixed")
AHR_roast_df <- as.matrix(rep(0,length(roast_cols)),nrow=1) %>% t() %>% as.data.frame()
colnames(AHR_roast_df) <- roast_cols
AHR_roast_df$NGenes <- AHR_roast$ngenes.in.set
AHR_roast_df$PropDown <- AHR_roast$p.value$Active.Prop[1] %>% round(.,4)
AHR_roast_df$PropUp <- AHR_roast$p.value$Active.Prop[2] %>% round(.,4)
AHR_roast_df$Direction <- ifelse(AHR_roast_df$PropDown > AHR_roast_df$PropUp, "Down", "Up")
AHR_roast_df$PValue <- ifelse(AHR_roast_df$Direction=="Down",
                              round(AHR_roast$p.value$P.Value[1], 4),
                              round(AHR_roast$p.value$P.Value[2],4))
AHR_roast_df$FDR <- round(AHR_roast$p.value$P.Value[3],4)
AHR_roast_df$PValue.Mixed <- AHR_roast_df$FDR.Mixed <- round(AHR_roast$p.value$P.Value[4],4)
AHR_roast_df 

write.table(AHR_roast_df, "../Results/Tables/Riaz_ICB_dataset_roast_res.txt", sep = "\t")

## calculate Immunophenoscore
ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}

## Assign colors
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
mapcolors<-function (x) {
  za<-NULL
  if (x>=3) {
    za=1000
  } else {
    if (x<=-3) {
      za=1
    } else {
      za=round(166.5*x+500.5,digits=0)
    }
  }
  return(my_palette[za])
}
my_palette2 <- colorRampPalette(c("black", "white"))(n = 1000)
mapbw<-function (x) {
  za2<-NULL
  if (x>=2) {
    za2=1000
  } else {
    if (x<=-2) {
      za2=1
    } else {
      za2=round(249.75*x+500.5,digits=0)
    }
  }
  return(my_palette2[za2])
}

## Estimate Immunophenoscore
IPS_function <- function(dge2,tcga_names){
  df <- data.frame(dge2)
  gene_expression <- df
  sample_names<-names(gene_expression)
  
  ## Read IPS genes and corresponding weights from tab-delimited text file "IPS_genes.txt"
  # For different 
  IPSG<-read.table("../Resources/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
  IPSG$GENE <- gsub("\\-","\\.",IPSG$GENE)
  IPSG$NAME <- gsub("\\-","\\.",IPSG$NAME)
  
  unique_ips_genes<-as.vector(unique(IPSG$NAME))
  
  IPS<-NULL
  MHC<-NULL
  CP<-NULL
  EC<-NULL
  SC<-NULL
  AZ<-NULL
  
  # Gene names in expression file
  GVEC<-row.names(gene_expression)
  # Genes names in IPS genes file
  VEC<-as.vector(IPSG$GENE)
  # Match IPS genes with genes in expression file
  ind<-which(is.na(match(VEC,GVEC)))
  # List genes missing or differently named
  MISSING_GENES<-VEC[ind]
  dat<-IPSG[ind,]
  if (length(MISSING_GENES)>0) {
    cat("differently named or missing genes: ",MISSING_GENES,"\n")
  }
  for (x in 1:length(ind)) {
    print(IPSG[ind,])
  }
  
  for (i in 1:length(sample_names)) {	
    GE<-gene_expression[[i]]
    mGE<-mean(GE)
    sGE<-sd(GE)
    Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
    W1<-IPSG$WEIGHT
    WEIGHT<-NULL
    MIG<-NULL
    k<-1
    for (gen in unique_ips_genes) {
      MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
      WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
      k<-k+1
    }
    WG<-MIG*WEIGHT
    MHC[i]<-mean(WG[1:10])
    CP[i]<-mean(WG[11:20])
    EC[i]<-mean(WG[21:24])
    SC[i]<-mean(WG[25:26])
    AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
    IPS[i]<-ipsmap(AZ[i])
  }
  
  df2 <- apply(df,2,function(x){(x-mean(x))/sd(x)}) %>% data.frame()
  il4i1 <- df2[rownames(df2) %in%unique(c("IL4I1", "IDO1", "TDO2")),]
  met_DF <- data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
  met_DF2 <- data.frame(met_DF,t(il4i1[,match(met_DF$SAMPLE,colnames(il4i1))]))
  met_DF2
}

dge_cpm <- cpm(dge$counts, prior.count = 1, log = T)

rownames(dge_cpm) <- dge$genes$genes

IPS_riaz <- IPS_function(dge_cpm, "raiz")

IPS_riaz2 <- data.frame(dge$samples,IPS_riaz)
IPS_riaz2$CP <- IPS_riaz2$CP*-1
IPS_riaz2$SC <- IPS_riaz2$SC*-1
IPS_riaz2$diff <- IPS_riaz2$EC-IPS_riaz2$SC
IPS_riaz2$PreOn <- factor(IPS_riaz2$PreOn, levels = c("Pre","On"))

## Indices of the paired patients
paired_pt_idcs <- IPS_riaz2$PatientID[which(duplicated(IPS_riaz2$PatientID)=="TRUE")]

df_plots_melt <- reshape::melt(IPS_riaz2)

df_paired <- df_plots_melt[df_plots_melt$PatientID %in% paired_pt_idcs,]
df_paired$Cohort[df_paired$Cohort=="NIV3-NAIVE"] <- "Ipilimumab_Naive"
df_paired$Cohort[df_paired$Cohort=="NIV3-PROG"] <- "Ipilimumab_Treated"

sel_cols <- c("IL4I1","IDO1","TDO2","MHC","CP","SC","EC","diff")

naive_df <- df_paired[df_paired$Cohort=="Ipilimumab_Naive",]%>%.[.$variable%in%sel_cols,]%>%
  .[.$PatientID!="Pt84",]%>% .[!is.na(.$Response),]

treated_df <- df_paired[df_paired$Cohort=="Ipilimumab_Treated",]%>%.[.$variable%in%sel_cols,]%>%
  .[.$PatientID!="Pt84",]%>% .[!is.na(.$Response),]

pdf("../Results/Figures/Riaz_ICB_dataset_naive_IL4I1_IDO1_CP_clinical_groups.pdf")
ggboxplot(naive_df%>% .[.$variable %in% c("IL4I1","IDO1","CP"),] ,
          x="PreOn", y="value", color = "PreOn",outlier.shape = NA,legend="top",
          add.params = list(color="black",size=2), size=2, title = "Naive", scales="free",
          add = "jitter", facet.by = c("variable","Response"), palette = c("blue","red"))+
  stat_compare_means(label = "p.format", paired = T, label.x = 1.3)
dev.off()

pdf("../Results/Figures/Riaz_ICB_dataset_treated_and_naive_IL4I1.pdf")
ggboxplot(df_paired %>% .[.$PatientID!="Pt84",]%>% .[!is.na(.$Response),] %>% .[.$variable %in% c("IL4I1"),], "PreOn","value", color="PreOn",
          facet.by = "Cohort", ncol=1, nrow=2, palette = c("blue","red"), add="jitter")+stat_compare_means(paired = T)
dev.off()

