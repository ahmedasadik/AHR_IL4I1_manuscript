# This script generates separate heatmaps portraying infiltrating cells across patients grouped by AHR activity

## Load libraries
library(gplots)
library(RColorBrewer)
library(purrr)

# The estimated GSVA scores of of the different cell types in all tumors
GSVA_voom_Immune_cells <- readRDS("./RDS/TCGA_GSVA_voom_28_cell_types_IPS.rds")

# Cluster Assignments
ccp_cluster_assignments <- readRDS("./RDS/ccp_cluster_assignments.rds")

# TCGA_clusters
TCGA_clusters <- paste("K",c(2,4,4,4,4,4,2,4,4,4,2,3,3,4,3,4,4,2,4,2,2,2,4,4,2,2,3,3,3,3,2,3), sep = "")

# tcga_names
tcga_names <- c("TCGA_ACC", "TCGA_BLCA", "TCGA_BRCA", "TCGA_CESC", "TCGA_COAD", "TCGA_CHOL", "TCGA_COAD", "TCGA_DLBC",
                "TCGA_ESCA", "TCGA_GBM", "TCGA_HNSC", "TCGA_KICH", "TCGA_KIRC", "TCGA_KIRP", "TCGA_LGG", "TCGA_LIHC",
                "TCGA_LUAD", "TCGA_LUSC", "TCGA_MESO", "TCGA_OV", "TCGA_PAAD", "TCGA_PCPG", "TCGA_PRAD", "TCGA_READ",
                "TCGA_SARC", "TCGA_SKCM", "TCGA_STAD", "TCGA_TGCT", "TCGA_THCA", "TCGA_THYM", "TCGA_UCEC", "TCGA_UCS", "TCGA_UVM")

# colors
cluster_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                    "#e377c2","#7f7f7f","#bcbd22","#17becf")

# heatmap colors
heatmap_colors <- brewer.pal(5,"RdBu")[5:1]

# Heatmap functions
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

## plot and save the heatmaps
walk(1:32,function(i,a1,a2,a3,a4){
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
  tiff(paste("./Figures/Heatmaps/GSVA_immune_cells", a3[i],a4[i],".tiff",sep = "_"),width = 12, height = 8, units = "in", res = 300)
  hm_RPPA_FUN(as.matrix(t(a1[[i]])), a2[[i]], ttl=paste(a3[i],a4[i],sep = "_"),
              sep_vec =sep_vec, col_p = colorRampPalette(heatmap_colors)(100))
  dev.off()
},a1=GSVA_voom_Immune_cells,a2=ccp_cluster_assignments,a3=tcga_names, a4=TCGA_clusters)
