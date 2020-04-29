# This script show how IDO1 and TDO2 strongly associate with AHR activation

## Load libraries
library(purrr)
library(edgeR)
library(corrplot)

## Functions
# Divide into groups based on median absolute deviation increments
group_no_z_FUN <- function(x,y){
  x_med <- median(x)
  x_mad <- mad(x)
  group <- vector("character", length = length(x))
  group[which(x >= (x_med + (y*x_mad)))] <- "high"
  group[which(x < (x_med - (y*x_mad)))] <- "low"
  group
}

# Perform gene set testing using roast (developed by Gordon Smyth)
GSA_roast_hi_lo_grps_FUN <- function(cnts,dge,goi,glist,sd_val){
  OV_TDO2_idx <- which(rownames(cnts)==goi)
  OV_TDO2_vals <- cnts[OV_TDO2_idx,]
  OV_med_0_group <- group_no_z_FUN(OV_TDO2_vals,sd_val)
  names(OV_med_0_group) <- colnames(cnts)
  if(sd_val==0){
    ct <- factor(OV_med_0_group, levels = c("low","high"))
    des_mat <- model.matrix(~ct) 
    index.vector <- colnames(dge$E) %in% glist
    AHR_roast <- roast(t(dge$E), index = index.vector, design = des_mat,
                       contrast = 2, set.statistic = "floormean")
  } else if(sd_val>0){
    OV_grps <- OV_med_0_group[-which(OV_med_0_group=="")]
    ct <- factor(OV_grps, levels = c("low","high"))
    des_mat <- model.matrix(~ct)
    index.vector <- colnames(dge$E) %in% glist
    AHR_roast <- roast(t(dge$E), index = index.vector, design = des_mat,
                       contrast = 2, set.statistic = "floormean")
  }
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

## Load Data
# TCGA_counts
TCGA_counts <- readRDS("./RDS/TCGA_counts.rds")

# TCGA_voom
TCGA_voom <- readRDS("./RDS/TCGA_DGE_voom_annot.rds")

# Read the AHR signature file
overlapping_genes <- read.delim("./Signature/overlapping_AHR_signature.txt", sep = "\t", stringsAsFactors = F)


####################################################################################
## AHR activation is associated with TDO2 and IDO1 expression in different tumors ##
####################################################################################
# Use the TCGA voom data set to avoid additional calculation steps, then run gene set testing

# TDO2 groups separated by the median into high and low. The comparison is high-low.
TDO2_GSA_roast_00_incr <- map2(TCGA_counts, TCGA_voom, GSA_roast_hi_lo_grps_FUN, "TDO2", overlapping_genes$Gene, 0) %>% do.call(rbind,.)
write.table(TDO2_GSA_roast_00_incr,"./Tables/TDO2_GSA_roast_00_incr.txt", sep = "\t")

# IDO1 groups separated by the median into high and low. The comparison is high-low.
IDO1_GSA_roast_00_incr <- map2(TCGA_counts, TCGA_voom, GSA_roast_hi_lo_grps_FUN, "IDO1", overlapping_genes$Gene, 0) %>% do.call(rbind,.)
write.table(IDO1_GSA_roast_00_incr,"./Tables/IDO1_GSA_roast_00_incr.txt", sep = "\t")

## Represent the up and down regulation as percentages of the proportional genes that are significantly different
GSA_corrplot_r <- data.frame(IDO1_dn=IDO1_GSA_roast_00_incr$PropDown*-1,
                             IDO1_up=IDO1_GSA_roast_00_incr$PropUp,
                             blank=rep(0,32),
                             TDO2_dn=TDO2_GSA_roast_00_incr$PropDown*-1,
                             TDO2_up=TDO2_GSA_roast_00_incr$PropUp)

GSA_corrplot_p <- data.frame(IDO1_dn=IDO1_GSA_roast_00_incr$FDR,
                             IDO1_up=IDO1_GSA_roast_00_incr$FDR,
                             blank=rep(1,32),
                             TDO2_dn=TDO2_GSA_roast_00_incr$FDR,
                             TDO2_up=TDO2_GSA_roast_00_incr$FDR)

rownames(GSA_corrplot_r) <- rownames(GSA_corrplot_p) <- gsub("TCGA_","",tcga_names)
colnames(GSA_corrplot_r) <- colnames(GSA_corrplot_p) <- c("IDO1_Dn", "IDO1_Up"," ", "TDO2_Dn", "TDO2_Up")

pdf("./Figures/GSA_IDO_TDO_median.pdf", width = 8, height = 12, pointsize = 12)
corrplot(corr = as.matrix(GSA_corrplot_r), p.mat = as.matrix(GSA_corrplot_p), insig = "blank", mar=c(0,0,1,0),
         cl.length = 3, col = colorRampPalette(colors = c("blue", "white", "red"))(11), addgrid.col = NA,
         tl.cex = 0.7, is.corr = F, method = "pie", cl.pos = "b", cl.lim = c(-1,1), outline = "black",
         tl.col = "black")
dev.off()
