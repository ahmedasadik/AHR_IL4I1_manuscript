# This script describes the validation of the separation between the clustered groups using PCA

# Load the CCP results
ccp_pos_dist <- readRDS("./RDS/ccp_pos_dist.rds")

## Load libraries
library(purrr)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(GGally)

## colors 
cluster_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                    "#e377c2","#7f7f7f","#bcbd22","#17becf")

# Read the overlapping AHR associated modules, these contain the Eigen-Genes
GT_GSVA_overlap_MEs <- readRDS("./RDS/GT_GSVA_overlap_MEs.rds")

# Results of concensus clustering
ccp_pos_dist <- readRDS("./RDS/ccp_pos_dist.rds")

## set seed
set.seed(0860)

## Running PCA
PCA_res <- purrr::map(GT_GSVA_overlap_MEs, PCA)

## PCA scree plots
PCA_Scree_plots <- purrr::map2(PCA_res, names(PCA_res), function(a,b){
  fviz_screeplot(a,main=b, addlabels = TRUE, linecolor="red", barcolor="grey")})

walk2(paste("./Figures/PCA_results/Plots_scree/",names(PCA_Scree_plots),".png", sep = ""), PCA_Scree_plots, function(a,b){
  png(a, res = 600, width = 12, height = 8, units = "in")
  print(b)
  dev.off()
})

### PCA paired plots. cluster colors are assigned 
PCA_pairs_plots <- purrr::map(1:2,function(idcs,a1, a2, a3){
  df <- as.data.frame(a1[[idcs]]$ind$coord)
  ccp_input <- a2[[idcs]]$ccp
  ccp_idcs <- cbind(2,3,4,5)
  ccp_idcs_mat <- apply(ccp_idcs,2, function(x){ccp_input[[x]]$consensusClass})
  paired_plots <- apply(ccp_idcs_mat, 2, function(ccp_i){
    p <- GGally::ggpairs(df, aes(col=as.factor(ccp_i)), title=a3[idcs])
    for(i in 1:p$nrow) {
      for(j in 1:p$ncol){
        p[i,j] <- p[i,j] + 
          scale_fill_manual(values=cluster_colors) +
          scale_color_manual(values=cluster_colors)+
          theme_bw()+theme(panel.grid = element_blank())
      }
    }
    p
  })
  names(paired_plots) <- c("K2", "K3", "K4", "K5")
  paired_plots
}, a1=PCA_res, a2=ccp_pos_dist, a3=names(ccp_pos_dist))
names(PCA_pairs_plots) <- names(PCA_res)[1:2]

## saving the paired plots to files
TCGA_PCA_tumor_idcs_for_saving <- rep(c(1:32), each=4)
PCA_pairs_plots_idcs_for_saving <- rep(c(1:4),32)
pairs_plots_to_save <- purrr::map2(TCGA_PCA_tumor_idcs_for_saving,PCA_pairs_plots_idcs_for_saving,
                                   function(a,b){p <- PCA_pairs_plots[[a]][b]})

to_save_names <- purrr::map2(TCGA_PCA_tumor_idcs_for_saving, PCA_pairs_plots_idcs_for_saving,
                             function(a,b){paste(names(PCA_res)[a], names(PCA_pairs_plots[[a]])[b], sep = "_")}) %>% unlist()

purrr::walk2(paste("./Figures/PCA_results/Plots_pairs/",to_save_names,".png", sep = ""), pairs_plots_to_save,
             function(a,b){
               png(a, res = 600, width = 12, height = 8, units = "in")
               print(b)
               dev.off()
             })

## PCA 3D plots
require(plotly)
PCA_3D_plots <- purrr::map(1:32,function(idcs,a1, a2, a3){
  df <- as.data.frame(a1[[idcs]]$ind$coord)
  ccp_input <- a2[[idcs]]$ccp
  ccp_idcs <- cbind(2,3,4,5)
  ccp_idcs_mat <- apply(ccp_idcs,2, function(x){
    ccp_input[[x]]$consensusClass})
  paired_plots <- apply(ccp_idcs_mat,2, function(ccp_i){
    p <- plot_ly(df, split=ccp_i, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, type = "scatter3d", mode = "lines")%>% layout(title = a3[[idcs]]) 
  }) 
  names(paired_plots) <- c("K2", "K3", "K4", "K5")
  paired_plots
}, a1=PCA_res, a2=ccp_pos_dist, a3=names(GT_GSVA_overlap_MEs))
names(PCA_3D_plots) <- names(GT_GSVA_overlap_MEs)

saveRDS(PCA_3D_plots, "./RDS/PCA_3D_plots.rds")
