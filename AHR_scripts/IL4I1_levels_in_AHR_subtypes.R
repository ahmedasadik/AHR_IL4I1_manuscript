# This scripts generates plots of IL4I1 expression in different tumors' AHR subtypes 

## Load libraries
library(purrr)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)

# TCGA_voom
TCGA_voom <- readRDS("/home/analyses/Projects_WF/WGCNA/TCGA_DGE_voom_annot.rds")

# Cluster Assignments
ccp_cluster_assignments <- readRDS("./RDS/ccp_cluster_assignments.rds")

# TCGA_clusters
TCGA_clusters <- paste("K",c(2,4,4,4,4,4,2,4,4,4,2,3,3,4,3,4,4,2,4,2,2,2,4,4,2,2,3,3,3,3,2,3), sep = "")

# colors
cluster_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                    "#e377c2","#7f7f7f","#bcbd22","#17becf")

# plot annotation function
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

# boxplot function
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

# Extract IL4I1 expression values
IL4I1_voom <- map(TCGA_voom, function(x)x$E[,which(colnames(x$E)=="IL4I1")])

# Create a dataframe with IL4I1 expression and cluster assignments
IL4I1_voom_boxplots <- map2(IL4I1_voom, ccp_cluster_assignments, function(a1,a2){
  df <- data.frame(IL4I1=a1, cluster=as.character(a2))
})

# Add the tumor type to the boxplot dataframe
IL4I1_voom_boxplots <- map2(IL4I1_voom_boxplots, names(TCGA_voom), function(a1,a2){
  df <- data.frame(tumor=a2, a1, stringsAsFactors = F)
})

# Generate plots
IL4I1_boxplots_list <- map2(IL4I1_voom_boxplots,names(TCGA_voom),dfs_to_plot_FUN)

# save plots
IL4I1_boxplots_arranged <- ggarrange(IL4I1_boxplots_list[[1]],IL4I1_boxplots_list[[2]],IL4I1_boxplots_list[[3]],
                                     IL4I1_boxplots_list[[4]],IL4I1_boxplots_list[[5]],IL4I1_boxplots_list[[6]],
                                     IL4I1_boxplots_list[[7]],IL4I1_boxplots_list[[8]],IL4I1_boxplots_list[[9]],
                                     IL4I1_boxplots_list[[10]],IL4I1_boxplots_list[[11]],IL4I1_boxplots_list[[12]],
                                     IL4I1_boxplots_list[[13]],IL4I1_boxplots_list[[14]],IL4I1_boxplots_list[[15]],
                                     IL4I1_boxplots_list[[16]],IL4I1_boxplots_list[[17]],IL4I1_boxplots_list[[18]],
                                     IL4I1_boxplots_list[[19]],IL4I1_boxplots_list[[20]],IL4I1_boxplots_list[[21]],
                                     IL4I1_boxplots_list[[22]],IL4I1_boxplots_list[[23]],IL4I1_boxplots_list[[24]],
                                     IL4I1_boxplots_list[[25]],IL4I1_boxplots_list[[26]],IL4I1_boxplots_list[[27]],
                                     IL4I1_boxplots_list[[28]],IL4I1_boxplots_list[[29]],IL4I1_boxplots_list[[30]],
                                     IL4I1_boxplots_list[[31]],IL4I1_boxplots_list[[32]], ncol = 4, nrow = 8)

pdf("./Figures/IL4I1_exprsn_all_AHR_groups_all_tumors.pdf",width = 18,height = 36,pointsize = 8)
IL4I1_boxplots_arranged
dev.off()
