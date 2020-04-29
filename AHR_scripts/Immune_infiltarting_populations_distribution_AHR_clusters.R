# This script describes creating a a series of heatmaps that show the different infiltrating immune
# populations in the AHR defined subtypes in the TCGA tumors

## Load libraries
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(RColorBrewer)

# The estimated GSVA scores of of the different cell types in all tumors
GSVA_voom_Immune_cells <- readRDS("./RDS/TCGA_GSVA_voom_28_cell_types_IPS.rds")

# Cluster Assignments
ccp_cluster_assignments <- readRDS("./RDS/ccp_cluster_assignments.rds")

# colors
cluster_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                    "#e377c2","#7f7f7f","#bcbd22","#17becf")

# heatmap colors
heatmap_colors <- brewer.pal(5,"RdBu")[5:1]

# Average the infiltration score per group across the tumors
GSVA_voom_Immune_cells_averaged_heatmaps <- map2(GSVA_voom_Immune_cells, ccp_cluster_assignments, function(a1,a2){
  df <- data.frame(cluster=a2, t(a1),  stringsAsFactors = F)
  split_levels <- levels(as.factor(df$cluster))
  split_dfs <- map(split_levels, function(b1){
    s_df <- df[df$cluster==b1,-1]
    colSums(s_df)
  })
  names(split_dfs) <- paste("G",split_levels, sep = "")
  res <- do.call(cbind,split_dfs)
})

# The defined AHR active clusters of all tumors that were determined from the roast run
heatmap_active_clusters <- c(2,4,1,3,2,3,2,2,3,1,2,1,2,1,1,1,3,2,1,2,2,1,3,4,2,1,3,1,3,1,2,3)

Averaged_GSVA_immune_heatmaps_all <- map(1:32, function(i,a1,a2,a3){
  annot_colors <- cluster_colors[1:ncol(a1[[i]])]
  names(annot_colors) <- colnames(a1[[i]])
  active_colors <- rep(cluster_colors[a3[i]], length(names(annot_colors)))
  names(active_colors) <- names(annot_colors)
  h1 <- HeatmapAnnotation(df = data.frame(cluster = names(annot_colors)),
                          col = list(cluster = annot_colors), show_legend = c(FALSE,FALSE))
  h2 <- HeatmapAnnotation(df = data.frame(Active = names(active_colors)),
                          col = list(Active = active_colors), show_legend = c(FALSE,FALSE))
  if(a2[i]=="ACC"){
    h3 <- Heatmap(a1[[i]], top_annotation = h1, bottom_annotation = h2, column_title = a2[i],
                  cluster_columns = FALSE, col = heatmap_colors, cluster_rows = FALSE,
                  row_names_side = "left",show_heatmap_legend = FALSE) 
  } else {
    h3 <- Heatmap(a1[[i]], show_row_names = FALSE,
                  top_annotation = h1,  bottom_annotation = h2, column_title = a2[i], cluster_columns = FALSE,
                  col = heatmap_colors, cluster_rows = FALSE,row_names_side = "left",show_heatmap_legend = FALSE)
  }
  h3
}, a1=GSVA_voom_Immune_cells_averaged_heatmaps,
a2=gsub("TCGA_","",names(GSVA_voom_Immune_cells_averaged_heatmaps)), a3=heatmap_active_clusters)

names(Averaged_GSVA_immune_heatmaps_all) <- names(TCGA_voom)

# Save a heatmap of all infiltration patterns across all the tumor clusters
pdf("./Figures/Immune_heatmap.pdf", width = 32, height = 8, pointsize = 18)
draw(Averaged_GSVA_immune_heatmaps_all[[1]]+Averaged_GSVA_immune_heatmaps_all[[2]]+Averaged_GSVA_immune_heatmaps_all[[3]]+
       Averaged_GSVA_immune_heatmaps_all[[4]]+Averaged_GSVA_immune_heatmaps_all[[5]]+Averaged_GSVA_immune_heatmaps_all[[6]]+
       Averaged_GSVA_immune_heatmaps_all[[7]]+Averaged_GSVA_immune_heatmaps_all[[8]]+Averaged_GSVA_immune_heatmaps_all[[9]]+
       Averaged_GSVA_immune_heatmaps_all[[10]]+Averaged_GSVA_immune_heatmaps_all[[11]]+Averaged_GSVA_immune_heatmaps_all[[12]]+
       Averaged_GSVA_immune_heatmaps_all[[13]]+Averaged_GSVA_immune_heatmaps_all[[14]]+Averaged_GSVA_immune_heatmaps_all[[15]]+
       Averaged_GSVA_immune_heatmaps_all[[16]]+Averaged_GSVA_immune_heatmaps_all[[17]]+Averaged_GSVA_immune_heatmaps_all[[18]]+
       Averaged_GSVA_immune_heatmaps_all[[19]]+Averaged_GSVA_immune_heatmaps_all[[20]]+Averaged_GSVA_immune_heatmaps_all[[21]]+
       Averaged_GSVA_immune_heatmaps_all[[22]]+Averaged_GSVA_immune_heatmaps_all[[23]]+Averaged_GSVA_immune_heatmaps_all[[24]]+
       Averaged_GSVA_immune_heatmaps_all[[25]]+Averaged_GSVA_immune_heatmaps_all[[26]]+Averaged_GSVA_immune_heatmaps_all[[27]]+
       Averaged_GSVA_immune_heatmaps_all[[28]]+Averaged_GSVA_immune_heatmaps_all[[29]]+Averaged_GSVA_immune_heatmaps_all[[30]]+
       Averaged_GSVA_immune_heatmaps_all[[31]]+Averaged_GSVA_immune_heatmaps_all[[32]], gap=unit(6, "mm"),show_annotation_legend = FALSE)
dev.off()
