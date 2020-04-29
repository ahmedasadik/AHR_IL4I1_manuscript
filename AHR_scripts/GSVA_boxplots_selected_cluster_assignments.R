# This script describes the generation of boxplots of the GSVA scores as per the selected cluster assignments

## Load libraries
library(purrr)
library(ggpubr)
library()

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

# TCGA_GSVA
TCGA_gsva <- readRDS("./RDS/TCGA_GSVA_scores_safely.rds")
# Extract the GSVA scores from the safely run output
TCGA_GSVA <- TCGA_gsva[tcga_names]
TCGA_GSVA <- map(TCGA_GSVA, function(x){x$result})

# Concensus clustering results
ccp_pos_dist <- readRDS("./RDS/ccp_pos_dist.rds")

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

## extract the CCP cluster assignments for each tumor after they were analyzed 
ccp_cluster_assignments <- map2(ccp_pos_dist,gsub("K","",TCGA_clusters)%>%as.numeric(), function(a1, a2){a1$ccp[[a2]]$consensusClass})
saveRDS(ccp_cluster_assignments,"./RDS/ccp_cluster_assignments.rds")

## Plot the boxplot representations
GSVA_clusters_plots <- map(1:32, function(i,a1,a2,a3){
  df_test <- data.frame(gsva_sc=t(a1[[i]]), clust=as.factor(a2[[i]]))
  comp_means_res <- compare_means(df_test, formula = gsva_sc~clust) %>%  as.data.frame()
  if(length(levels(df_test$clust))==2){
    comp_y_position <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
  } else if(length(levels(df_test$clust))==3){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
    comp_y_position <- seq(from = val_from,by = 0.2,length.out = 3)
  } else if(length(levels(df_test$clust))==4){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
    comp_y_position <- seq(from = val_from,by = 0.2,length.out = 6)
  }
  comp_means_res <-comp_means_res %>% mutate(.,y.position=comp_y_position) %>%
    mutate(sig=annt_fun(.$p.adj,length(.$p.adj)))
  p <- ggboxplot(df_test, x="clust", y="gsva_sc", color = "clust",
                 legend="right",palette = cluster_colors,ticks=FALSE,
                 scales="free_y",x.text.angle = 45,
                 title = a3[[i]])+
    stat_pvalue_manual(label = "p.adj",comp_means_res)+
    ggbeeswarm::geom_quasirandom(aes(color=clust))
  p
}, a1=TCGA_GSVA, a2=ccp_cluster_assignments, a3=tcga_names)

names(GSVA_clusters_plots) <- tcga_names

## Save
walk2(GSVA_clusters_plots, paste("./Figures/GSVA_boxplots_ccp/GSVA_clust_", tcga_names, ".png", sep = ""),
      function(a1,a2){
        png(a2, res = 600, width = 12, height = 8, units = "in")
        print(a1)
        dev.off()
      })
