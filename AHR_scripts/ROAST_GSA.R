# This script describes running roast to evaluate the up or down regulation of AHR activity among different groups

## Load libraries
library(purrr)
library(edgeR)

## Read the AHR signature file
overlapping_genes <- read.delim("./Signature/overlapping_AHR_signature.txt", sep = "\t", stringsAsFactors = F)

# TCGA_voom
TCGA_voom <- readRDS("./RDS/TCGA_DGE_voom_annots.rds")

# Concensus clustering results
ccp_pos_dist <- readRDS("./RDS/ccp_pos_dist.rds")

# selected TCGA_clusters comparisons toptables
TCGA_DGE_ccp_clusters <- readRDS("./RDS/TCGA_DGE_ccp_clusters.rds")

# The selected number of clusters of each TCGA tumor (this was analyzed manually as previously described)
TCGA_clusters <- paste("K",c(2,4,4,4,4,4,2,4,4,4,2,3,3,4,3,4,4,2,4,2,2,2,4,4,2,2,3,3,3,3,2,3), sep = "")

# tcga_names
tcga_names <- c("TCGA_ACC", "TCGA_BLCA", "TCGA_BRCA", "TCGA_CESC", "TCGA_COAD", "TCGA_CHOL", "TCGA_COAD", "TCGA_DLBC",
                "TCGA_ESCA", "TCGA_GBM", "TCGA_HNSC", "TCGA_KICH", "TCGA_KIRC", "TCGA_KIRP", "TCGA_LGG", "TCGA_LIHC",
                "TCGA_LUAD", "TCGA_LUSC", "TCGA_MESO", "TCGA_OV", "TCGA_PAAD", "TCGA_PCPG", "TCGA_PRAD", "TCGA_READ",
                "TCGA_SARC", "TCGA_SKCM", "TCGA_STAD", "TCGA_TGCT", "TCGA_THCA", "TCGA_THYM", "TCGA_UCEC", "TCGA_UCS", "TCGA_UVM")

## Create design matrices
TCGA_des_mats <- map2(ccp_pos_dist, TCGA_clusters, function(a1,a2){
  idx <- gsub("K","",a2) %>% as.numeric()
  clss <- a1$ccp[[idx]]$consensusClass
  des<-model.matrix(~0+as.factor(clss))
  colnames(des) <- paste("cls",levels(as.factor(clss)),sep="_")
  des})

## Create contrast matrices
TCGA_cont_mats <- map(TCGA_des_mats, function(des_mat){
  if(length(colnames(des_mat))==2){
    cont_mat <- makeContrasts(cls_1-cls_2, levels = des_mat)
  } else if(length(colnames(des_mat))==3){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_2-cls_3, levels = des_mat)
  } else if(length(colnames(des_mat))==4){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_1-cls_4,
                              cls_2-cls_3, cls_2-cls_4, cls_3-cls_4,levels = des_mat)
  }
})

## Run roast for the selected clusters ####
TCGA_roast_FUN <- function(idx, v_eset, glist, des_mat, cont_mat){
  index.vector <- colnames(v_eset$E) %in% glist
  AHR_roast <- roast(t(v_eset$E), index = index.vector, design = des_mat,
                     contrast = cont_mat[,idx], set.statistic = "floormean")
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

TCGA_roast_res <- map(1:32, function(i, v_eset, glist, des_mat, cont_mat){
  if(length(colnames(cont_mat[[i]])) == 1){
    cnts <- TCGA_roast_FUN(idx = 1, v_eset = v_eset[[i]], glist = glist, des_mat = des_mat[[i]], cont_mat[[i]])  
  } else {
    cnts <- map(1:length(colnames(cont_mat[[i]])),TCGA_roast_FUN,
                v_eset = v_eset[[i]], glist = glist, des_mat = des_mat[[i]], cont_mat[[i]]) 
  }
  cnts
}, v_eset=TCGA_voom, glist=overlapping_genes$Gene, des_mat=TCGA_des_mats, cont_mat=TCGA_cont_mats)

names(TCGA_roast_res) <- tcga_names

cnts_names <- map2(TCGA_DGE_ccp_clusters, names(TCGA_DGE_ccp_clusters), function(a1,a2){
  cnts_names <- names(a1)
  paste(a2, cnts_names, sep = "_")
}) %>% unlist()

TCGA_roast_res_df <- do.call(bind_rows, TCGA_roast_res) %>% data.frame(Contrasts=cnts_names,.,stringsAsFactors=FALSE)
TCGA_roast_res_df$pct_diff <- (abs(TCGA_roast_res_df$PropDown-TCGA_roast_res_df$PropUp)/(TCGA_roast_res_df$PropDown+TCGA_roast_res_df$PropUp))*100

write.csv(TCGA_roast_res_df, "./Tables/TCGA_roast_results.csv")

## After roast, order the groups based on AHR activation for the different tumor types ####
TCGA_DGE_ccp_clusters_AHR_activation_order <- map(gsub("TCGA_","",tcga_names), function(a1){
  df <- rep(NA,4) %>% as.data.frame() %>% t() %>% as.data.frame()
  colnames(df) <- paste("C",1:4, sep = "")
  rownames(df) <- a1
  idcs <- grep(a1, TCGA_roast_res_df$Contrasts)
  if(length(idcs)==1){
    if(TCGA_roast_res_df$Direction[idcs]=="Up"){
      df <- c(1,2,NA,NA)
    } else {
      df <- c(2,1,NA,NA)
    }
  } else if(length(idcs)==3){
    if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up"))){
      df <- c(1,2,3,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Down"))){
      df <- c(1,3,2,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Up"))){
      df <- c(2,1,3,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Up"))){
      df <- c(2,3,1,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Down"))){
      df <- c(3,1,2,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down"))){
      df <- c(3,2,1,NA)
    }
  } else if(length(idcs)==6){
    if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Up","Up","Up"))){
      df <- c(1,2,3,4)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Up","Up","Down"))){
      df <- c(1,2,4,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Down","Up","Up"))){
      df <- c(1,3,2,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Down","Down","Up"))){
      df <- c(1,3,4,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Up","Down","Down"))){
      df <- c(1,4,2,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Down","Down","Down"))){
      df <- c(1,4,3,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Up","Up","Up","Up"))){
      df <- c(2,1,3,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Up","Up","Up","Down"))){
      df <- c(2,1,4,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Up","Up","Up","Up"))){
      df <- c(2,3,1,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Up","Up","Up"))){
      df <- c(2,3,4,1)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Down","Up","Up","Down"))){
      df <- c(2,4,1,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Up","Up","Down"))){
      df <- c(2,4,3,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Up","Down","Up","Up"))){
      df <- c(3,1,2,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Up","Down","Down","Up"))){
      df <- c(3,1,4,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Up","Down","Up","Up"))){
      df <- c(3,2,1,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Down","Up","Up"))){
      df <- c(3,2,4,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Down","Down","Down","Up"))){
      df <- c(3,4,1,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Down","Down","Up"))){
      df <- c(3,4,2,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Down","Up","Down","Down"))){
      df <- c(4,1,2,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Down","Down","Down","Down"))){
      df <- c(4,1,3,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Down","Up","Down","Down"))){
      df <- c(4,2,1,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Up","Down","Down"))){
      df <- c(4,2,3,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Down","Down","Down","Down"))){
      df <- c(4,3,1,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Down","Down","Down"))){
      df <- c(4,3,2,1)
    }
    df
  }
}) %>% do.call(rbind,.)

TCGA_DGE_ccp_clusters_AHR_activation_order_df <- data.frame(Tumor=tcga_names, TCGA_DGE_ccp_clusters_AHR_activation_order)

write.table(TCGA_DGE_ccp_clusters_AHR_activation_order_df, "./Tables/TCGA_DGE_ccp_clusters_AHR_activation_order_df.txt")

## TCGA_barcode_plots for the selected clusters that were roasted
walk(1:32, function(i, tt_list, tcga_names, glist){
  if(length(names(tt_list[[i]])) == 1){
    index.vector2 <- tt_list[[i]][[1]]$ID %in% glist
    png(paste("./Figures/barcodeplots_", tcga_names[i],"_cnt1.png", sep=""),width = 12, height = 8, units = "in", res = 600)
    barcodeplot(tt_list[[i]][[1]]$t, index = index.vector2, main=paste(tcga_names[i],"cnt1", sep = "_"),
                gene.weights = tt_list[[i]][[1]][index.vector2,"logFC"])
    dev.off()
  } else {
    walk(1:length(names(tt_list[[i]])),function(idx){
      index.vector2 <- tt_list[[i]][[idx]]$ID %in% glist
      png(paste("./Figures/barcodeplots_", tcga_names[i],"_",names(tt_list[[i]][idx]),".png", sep=""),width = 12, height = 8, units = "in", res = 600)
      barcodeplot(tt_list[[i]][[idx]]$t, index = index.vector2, main=paste(tcga_names[i], names(tt_list[[i]][idx]), sep = "_"),
                  gene.weights = tt_list[[i]][[idx]][index.vector2,"logFC"])
      dev.off()
    })
  }
}, tt_list=TCGA_DGE_ccp_clusters, tcga_names=tcga_names, glist=overlapping_genes$Gene)

