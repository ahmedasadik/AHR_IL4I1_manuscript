# This script describes the generation of heatmaps for the different clustered groups
# using all modules and the genes of the modules associating only positively with AHR

## Load libraries
library(purrr)
library(gplots)
library(RColorBrewer)

# Heatmap function
hm_RPPA_FUN <- function(tt_df, clust_grp, ttl, col_p, sep_vec,...){
  col_clust <- clust_grp
  ordered_patients <- tt_df[order(col_clust),]
  to_plot <- t(scale(as.matrix(ordered_patients)))
  
  heatmap.2(to_plot,col=col_p, trace = "none",Colv = NULL,main = ttl, colsep = sep_vec,
            scale = "row", cexRow =0.8, dendrogram = "none", key.title = NA,margins = c(8, 8),
            density.info = c("density"), labRow = rownames(to_plot), labCol = FALSE, 
            ColSideColors = cluster_colors[col_clust[order(col_clust)]])
}

# TCGA_voom
TCGA_voom <- readRDS("./RDS/TCGA_DGE_voom_annots.rds")

# Read the overlapping AHR associated modules, these contain the Eigen-Genes
GT_GSVA_overlap_MEs <- readRDS("./RDS/GT_GSVA_overlap_MEs.rds")

# Pearson correlation of modules with AHR activation
TCGA_MEs_GSVA_cors <- readRDS("./RDS/TCGA_MEs_GSVA_cors.rds")

# Genes of WGCNA modules
TCGA_modules_genes <- readRDS("./RDS/TCGA_modules_genes.rds")

# Consensus cluster assignments
ccp_cluster_assignments <- readRDS("./RDS/ccp_cluster_assignments.rds")

# The selected number of clusters of each TCGA tumor (this was analyzed manually as previously described)
TCGA_clusters <- paste("K",c(2,4,4,4,4,4,2,4,4,4,2,3,3,4,3,4,4,2,4,2,2,2,4,4,2,2,3,3,3,3,2,3), sep = "")

# tcga_names
tcga_names <- c("TCGA_ACC", "TCGA_BLCA", "TCGA_BRCA", "TCGA_CESC", "TCGA_COAD", "TCGA_CHOL", "TCGA_COAD", "TCGA_DLBC",
                "TCGA_ESCA", "TCGA_GBM", "TCGA_HNSC", "TCGA_KICH", "TCGA_KIRC", "TCGA_KIRP", "TCGA_LGG", "TCGA_LIHC",
                "TCGA_LUAD", "TCGA_LUSC", "TCGA_MESO", "TCGA_OV", "TCGA_PAAD", "TCGA_PCPG", "TCGA_PRAD", "TCGA_READ",
                "TCGA_SARC", "TCGA_SKCM", "TCGA_STAD", "TCGA_TGCT", "TCGA_THCA", "TCGA_THYM", "TCGA_UCEC", "TCGA_UCS", "TCGA_UVM")

## colors
cluster_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                    "#e377c2","#7f7f7f","#bcbd22","#17becf")

## Heatmaps of all modules associating with AHR activation
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
  png(paste("./Figures/Heatmaps/MEs/MEs",a3[i],a4[i],".png",sep = "_"),width = 12, height = 8, units = "in", res = 600)
  hm_RPPA_FUN(as.matrix((a1[[i]])), a2[[i]], ttl=paste(a3[i],a4[i],sep = "_"), col_p=greenred(100),
              sep_vec =sep_vec)
  dev.off()
},a1=GT_GSVA_overlap_MEs, a2=ccp_cluster_assignments, a3=tcga_names, a4=TCGA_clusters)


## Heatmaps of genes of modules only positively associating with AHR activation ##

## From the GlobalTest results retain only the modules that have a positive significant
## association with AHR activation
TCGA_MEs_GSVA_GTs_pos_only <- map(GT_GSVA_overlap_MEs, function(x){
  x[grep("pos",x$direction),]
})

## Overlap the retained modules with the GSVA correlations 
TCGA_MEs_GSVA_GT_overlap_pos_only <- map2(TCGA_MEs_GSVA_GTs_pos_only, TCGA_MEs_GSVA_cors, function(x,y){
  GT_MEs <- gsub("ME","",rownames(x))
  Cor_MEs <- y$modules
  GT_MEs[match(Cor_MEs, GT_MEs)] %>% .[!is.na(.)]
})

## Extract the genes of the corresponding pos_only modules
TCGA_modules_genes_pos_only <- map2(TCGA_modules_genes, TCGA_MEs_GSVA_GT_overlap_pos_only, function(a1,a2){
  a1[a1$module %in% a2,]
})

TCGA_modules_genes_pos_only_cpm_mats <- map2(TCGA_voom, TCGA_modules_genes_pos_only,
                                             function(x,b){(x$E%>%t())[rownames(x$E%>%t()) %in% b$genes,]})

saveRDS(TCGA_modules_genes_pos_only_cpm_mats,"./RDS/TCGA_modules_genes_pos_only_cpm_mats.rds")

# heatmaps
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
  png(paste("./Figures/Heatmaps/ME_genes/ME_genes",a3[i],a4[i],".png",sep = "_"),width = 12, height = 8, units = "in", res = 600)
  hm_RPPA_FUN(as.matrix(t(a1[[i]])), a2[[i]], ttl=paste(a3[i],a4[i],sep = "_"), col_p=greenred(100),
              sep_vec =sep_vec, labRow=NA)
  dev.off()
},a1=TCGA_modules_genes_pos_only_cpm_mats, a2=ccp_cluster_assignments, a3=tcga_names, a4=TCGA_clusters)
