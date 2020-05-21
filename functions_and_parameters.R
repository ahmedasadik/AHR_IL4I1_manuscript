#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes the functions and parameters needed
################################################

#################################################
## VARIABLES
#################################################

my_seed <- 0860
no_cores <- parallel::detectCores() -1

# Selected enzymes of the Tryptophan degradation pathway
trp_enz_sel <- c("IL4I1", "IDO1", "IDO2", "TDO2", "DDC", "TPH1", "TPH2")

# tcga_project_ids
project_ids <- c("TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC",
                "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LGG",
                "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", 
                "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA",
                "TCGA-THYM", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM")
# tcga_names
tcga_names <- c("TCGA_ACC", "TCGA_BLCA", "TCGA_BRCA", "TCGA_CESC", "TCGA_CHOL", "TCGA_COAD", "TCGA_DLBC",
                "TCGA_ESCA", "TCGA_GBM", "TCGA_HNSC", "TCGA_KICH", "TCGA_KIRC", "TCGA_KIRP", "TCGA_LGG",
                "TCGA_LIHC", "TCGA_LUAD", "TCGA_LUSC", "TCGA_MESO", "TCGA_OV", "TCGA_PAAD", "TCGA_PCPG", 
                "TCGA_PRAD", "TCGA_READ", "TCGA_SARC", "TCGA_SKCM", "TCGA_STAD", "TCGA_TGCT", "TCGA_THCA",
                "TCGA_THYM", "TCGA_UCEC", "TCGA_UCS", "TCGA_UVM")

TCGA_names <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
                "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC",
                "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

#################################################
## FUNCTIONS
#################################################

# Annotation function used with HGNC input
annot_query_FUN <- function(gS, sym, eid, syn, psym=NA){
  if (length(grep(paste("^(", gS, ")$", sep = ""), sym)) > 0){
    idx <- grep(paste("^(", gS, ")$", sep = ""), sym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("^(", gS, ")$", sep = ""), psym)) > 0){
    idx <- grep(paste("^(", gS, ")$", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("^(", gS, ")\\|", sep = ""), psym)) > 0){
    idx <- grep(paste("\\|(", gS, ")\\|", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("\\|(", gS, ")\\|", sep = ""), psym)) > 0){
    idx <- grep(paste("\\|(", gS, ")\\|", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("\\|(", gS, ")$", sep = ""), psym)) > 0){
    idx <- grep(paste("\\|(", gS, ")$", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (length(grep(paste("^(", gS, ")$", sep = ""), syn)) > 0){
    idx <- grep(paste("^(", gS, ")$", sep = ""), syn)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else {
    idx <- grep(paste("^(", gS, ")\\|", sep = ""), syn)
    idx2 <- grep(paste("\\|(", gS, ")\\|", sep = ""), syn)
    idx3 <- grep(paste("\\|(", gS, ")$", sep = ""), syn)
    if (length(idx)>0){
      return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
    } else if (length(idx2)>0) {
      return(c(paste(sym[idx2], collapse = ";"), paste(eid[idx2], collapse = ";")))
    } else {
      return(c(paste(sym[idx3], collapse = ";"), paste(eid[idx3], collapse = ";")))
    }
  }
}

# Replace blank cells with NA 
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

# Functions to download and process the TCGA data and save it as an RData file
counts_download_FUN <- function(pid){
  query <- TCGAbiolinks::GDCquery(project = pid,
                                  data.category = "Transcriptome Profiling",
                                  experimental.strategy = "RNA-Seq",
                                  legacy = FALSE,
                                  workflow.type = "HTSeq - Counts")
  TCGAbiolinks::GDCdownload(query, method = "client")
  TCGAbiolinks::GDCprepare(query, save = T,
                           save.filename=paste0("./RDats/", query$project, gsub(" ", "_", query$data.category),gsub(" ","_",date()),".RData"))
}

fpkm_download_FUN <- function(pid){
  query <- TCGAbiolinks::GDCquery(project = pid,
                                  data.category = "Transcriptome Profiling",
                                  experimental.strategy = "RNA-Seq",
                                  legacy = FALSE,
                                  workflow.type = "HTSeq - FPKM")
  TCGAbiolinks::GDCdownload(query)
  TCGAbiolinks::GDCprepare(query, save = T,
                           save.filename=paste0("./RDats/", query$project, gsub(" ", "_", query$data.category),gsub(" ","_",date()),".RData"))
}
  
# Function that converts FPKMs to TPMs
fpkmToTpm <- function(fpkm){(fpkm / sum(fpkm)) * 1e6}

## Plot annotation function
annt_fun <- function(sig_diff_s, l){
  lbl <- vector("character", length = l)
  for (i in seq_along(sig_diff_s)){
    if(sig_diff_s[i] > 0.05){
      lbl[i] <- "ns"
    } else if (sig_diff_s[i] < 0.05 & sig_diff_s[i] > 0.01){
      lbl[i] <- "*"
    } else if (sig_diff_s[i] <= 0.01 & sig_diff_s[i] > 0.001){
      lbl[i] <- "**"
    } else if (sig_diff_s[i] <= 0.001 & sig_diff_s[i] > 0.0001){
      lbl[i] <- "***"
    } else if (sig_diff_s[i] <= 0.0001){
      lbl[i] <- "****"
    }
  }
  lbl
}

## Plot WGCNA circos plots
WGCNA_Circ_plot_FUN <- function(i, g_mats,  mod_w_enz, ttl, save_path, enz_vec){
  ## Color function
  col_fun = colorRamp2(c(-8, 0, 8), c("blue", "white", "red"))
  ## matrix for the run
  mat2=g_mats[[i]]
  enz_mat <- matrix(0,ncol=ncol(mat2), nrow = length(enz_vec)*50) %>% as.data.frame()
  enz_mat[,1] <- rep("trp_enz", length(enz_vec)*50)
  enz_mat[,2] <- rep(enz_vec, each=50)
  colnames(enz_mat) <- colnames(mat2)
  mat2 <- rbind(mat2, enz_mat) %>% data.frame(., stringsAsFactors=F)
  mat3 <- mat2[,-c(1:2)] %>% t()
  factors2 <- mat2$module
  to_subset <- as.character(levels(factors2))
  mat_list2.0 <- map(to_subset, function(x){mat3[,factors2==x]})
  names(mat_list2.0) <- to_subset
  mod_enz_pr <- mod_w_enz[[i]] %>% .[!is.na(.$aaa_enzymes),]
  pdf(save_path[i], width = 4, height = 4, pointsize = 16)
  ## plot outline
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 4)
  x_limits <- cbind(rep(0,length(to_subset)), table(factors2))
  circos.initialize(to_subset, xlim = x_limits)
  circos.track(ylim = c(0, 0.5), bg.border = NA, track.height=0.25,
               panel.fun = function(x, y) {
                 x_p <- get.cell.meta.data("xlim")
                 y_p <- get.cell.meta.data("ylim")
                 circos.rect(x_p[1], y_p[1],x_p[2], y_p[2],
                             col = ifelse(CELL_META$sector.index %in% colors(),CELL_META$sector.index, "white"),
                             border = ifelse(CELL_META$sector.index %in% colors(),"black", "white"))
               })
  ## add enzyme names
  map(enz_vec, function(z){
    circos.text(25, 0.35, labels=z, sector.index = z, cex=0.9,
                facing = "clockwise", niceFacing = TRUE, col = "black")
  })
  ## plot links
  mod_enz_pr <- mod_enz_pr[mod_enz_pr$genes%in%enz_vec,]
  if(dim(mod_enz_pr)[1]!=0){
    map2(mod_enz_pr$genes,mod_enz_pr$module, function(x,y){
      l_inv <- x_limits[rownames(x_limits)==y,2]/2
      circos.link(x, c(0,50), y, c(l_inv-10, l_inv+10), col = y, border = "black") 
    }) 
  }
  
  ## add title
  title(ttl[i])
  circos.clear()
  dev.off()
}

# Divide into groups based on median absolute deviation increments
group_no_z_FUN <- function(x,y){
  x_med <- median(x)
  x_mad <- mad(x)
  group <- vector("character", length = length(x))
  group[which(x >= (x_med + (y*x_mad)))] <- "high"
  group[which(x < (x_med - (y*x_mad)))] <- "low"
  group
}

# Perform gene set testing using roast (PMID: 20610611)
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

# Plotting functions
# plot annotation function
dfs_to_plot_FUN <- function(df,ttl){
  df_melt <- reshape::melt(df)
  comp_means_res <- compare_means(formula = value~cluster, data=df_melt) %>%  as.data.frame()
  if(length(levels(df_melt$cluster))==2){
    comp_y_position <- max(df_melt$value,na.rm = TRUE)+0.1
  } else if(length(levels(df_melt$cluster))==3){
    val_from <- max(df_melt$value,na.rm = TRUE)+0.1
    comp_y_position <- seq(from = val_from,by = 0.6,length.out = 3)
  } else if(length(levels(df_melt$cluster))==4){
    val_from <- max(df_melt$value,na.rm = TRUE)+0.1
    comp_y_position <- seq(from = val_from,by = 0.6,length.out = 6)
  }
  comp_means_res <-comp_means_res %>% mutate(.,y.position=comp_y_position) %>%
    mutate(sig=annt_fun(.$p.adj,length(.$p.adj)))
  p <- ggboxplot(df_melt, x="cluster", y="value", color = "cluster", facet.by = "variable",
                 legend="right",palette = cluster_colors, ticks=FALSE,ylab = ttl,
                 scales="free_y",x.text.angle = 45
  )+stat_pvalue_manual(label = "p.adj", comp_means_res)+ggbeeswarm::geom_quasirandom(aes(color=cluster))
  p
}

hm_RPPA_FUN <- function(tt_df, clust_grp, ttl, col_p, sep_vec,...){
  col_clust <- clust_grp
  ordered_patients <- tt_df[order(col_clust),]
  to_plot <- t(scale(as.matrix(ordered_patients)))
  
  heatmap.2(to_plot,col=col_p, trace = "none", Colv = NULL,Rowv = FALSE,main = ttl, colsep = sep_vec,sepcolor = "black",
            scale = "row", cexRow =0.8, dendrogram = "none", key.title = NA,margins = c(8, 8),
            density.info = c("density"), labRow = rownames(to_plot), labCol = FALSE,
            ColSideColors = cluster_colors[col_clust[order(col_clust)]])
}

HM_plot_FUN <- function(i,a1,a2,a3,a4){
  col_clust <- a2[[i]]
  grp_nums <- table(col_clust)
  sep_vec <- c()
  if(dim(grp_nums)==2){
    sep_vec <- grp_nums[1]
  } else if(dim(grp_nums)==3){
    sep_vec <- c(grp_nums[1],grp_nums[1]+grp_nums[2])
  } else if(dim(grp_nums)==4){
    sep_vec <- c(grp_nums[1],grp_nums[1]+grp_nums[2],grp_nums[1]+grp_nums[2]+grp_nums[3])
  }
  hm_RPPA_FUN(as.matrix(t(a1[[i]])), a2[[i]], ttl=paste(a3[i],a4[i],sep = "_"),
              col_p=colorRampPalette(heatmap_colors)(100),
              sep_vec =sep_vec)
}

avoid_overlap <- function(x){
  ind <- seq_along(x) %% 2==0
  x[ind] <- gsub("\\s", " ", format(x[ind], width=max(nchar(x[ind]))+8))
  x[ind] <- paste0(x[ind], "         ")
  x
}