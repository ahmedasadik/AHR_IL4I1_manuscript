# This script describes the generation of the WGCNA Circos-plot of GBM with the selected enzymes of the Trp degradation pathway

## Load libraries
library(purrr)
library(ComplexHeatmap)
library(circlize)

## WGCNA circos plot function
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
  tiff(save_path[i], width = 4, height = 4, units = "in", res = 300)
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
                             border = NA)
               })
  ## add enzyme names
  map(enz_vec, function(z){
    circos.text(25, 0.2, labels=z, sector.index = z, cex=0.8,
                facing = "clockwise", niceFacing = TRUE, col = "black")
  })
  ## plot links
  mod_enz_pr <- mod_enz_pr[mod_enz_pr$genes%in%enz_vec,]
  if(dim(mod_enz_pr)[1]!=0){
    map2(mod_enz_pr$genes,mod_enz_pr$module, function(x,y){
      l_inv <- x_limits[rownames(x_limits)==y,2]/2
      circos.link(x, c(0,50), y, c(l_inv-10, l_inv+10), col = y) 
    }) 
  }
  
  ## add title
  title(ttl[i])
  circos.clear()
  dev.off()
}

# TCGA_voom
TCGA_voom <- readRDS("./RDS/TCGA_DGE_voom_annots.rds")

# Genes of the WGCNA modules
TCGA_modules_genes <- readRDS("./RDS/TCGA_modules_genes.rds")

# Positive only correlating modules
TCGA_MEs_GSVA_cors_pos_only_overlap_GT <- readRDS("./RDS/TCGA_MEs_GSVA_cors_pos_only_overlap_GT.rds")

# Lists defining the presence of enzymes in the positive only associating modules
modules_with_enzymes_pos <- readRDS("./RDS/modules_with_enzymes_pos.rds")

# tcga_names
tcga_names <- c("TCGA_ACC", "TCGA_BLCA", "TCGA_BRCA", "TCGA_CESC", "TCGA_COAD", "TCGA_CHOL", "TCGA_COAD", "TCGA_DLBC",
                "TCGA_ESCA", "TCGA_GBM", "TCGA_HNSC", "TCGA_KICH", "TCGA_KIRC", "TCGA_KIRP", "TCGA_LGG", "TCGA_LIHC",
                "TCGA_LUAD", "TCGA_LUSC", "TCGA_MESO", "TCGA_OV", "TCGA_PAAD", "TCGA_PCPG", "TCGA_PRAD", "TCGA_READ",
                "TCGA_SARC", "TCGA_SKCM", "TCGA_STAD", "TCGA_TGCT", "TCGA_THCA", "TCGA_THYM", "TCGA_UCEC", "TCGA_UCS", "TCGA_UVM")

# Selected enzymes of the Tryptophan degradation pathway
trp_enz_sel <- c("IL4I1","IDO1","TDO2","KMO","KYNU","MAOB","MAOA","CAT","TPH1","TPH2","DDC","KYAT3","KYAT1","AADAT","AFMID")

save_paths <- paste("./Figures/WGCNA_circos_selected", tcga_names, ".tiff", sep = "")

## Matrices with the expression values of the genes in the different modules that are positively associated with AHR
TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules <- map(1:length(TCGA_modules_genes),function(i,mod_genes, gsva_corrs, voom_genes){
  res_cmg <- mod_genes[[i]][mod_genes[[i]]$module %in% gsva_corrs[[i]]$modules,]
  res_cmg_nums <- as.data.frame(table(res_cmg$module))
  res_cmg_nums <- res_cmg_nums[order(res_cmg_nums$Freq, decreasing = T),]
  res_cmg$module <- factor(res_cmg$module, levels = res_cmg_nums$Var1)
  res_cmg <- res_cmg[order(res_cmg$module),]
  res_cmg_E <- t(voom_genes[[i]]$E[,match(res_cmg$genes,colnames(voom_genes[[i]]$E))])
  res_cmg_all <- data.frame(res_cmg, res_cmg_E)
},
mod_genes=TCGA_modules_genes, gsva_corrs=TCGA_MEs_GSVA_cors_pos_only_overlap_GT, voom_genes=TCGA_voom)

names(TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules) <- tcga_names

WGCNA_Circ_plot_FUN(9, g_mats=TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules,
                    mod_w_enz=modules_with_enzymes_pos, ttl=gsub("TCGA_","",tcga_names),
                    save_path=save_paths, enz_vec=trp_enz_sel)
