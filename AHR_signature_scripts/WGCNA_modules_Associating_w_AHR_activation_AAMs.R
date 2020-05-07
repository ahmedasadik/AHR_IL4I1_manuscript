#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes how modules associating with AHR activation are defined.
## This is achieved by testing the association of AHR activation (GSVA scores) with the WGCNA modules.
## To overcome bias, we will try two different methods and take the overlap between them.
## Following this, circos plots are generated showing the presence of Trp catabolizing enzymes in the WGCNA modules with positive associations with AHR activation (AAMs)
#################################################

## Load libraries
library(purrr)
library(globaltest)
library(circlize)
library(ComplexHeatmap)

## Source the functions and parameters files
source("../functions_and_parameters.R")

## Data
# TCGA_modules
TCGA_toms <- readRDS("../Zenodo_download/TCGA_TOMS_bicor.rds")

# TCGA_modules_genes
TCGA_modules_genes <- readRDS("../Zenodo_download/TCGA_modules_genes.rds")

# TCGA_voom
TCGA_voom <- readRDS("../Zenodo_download/TCGA_DGE_voom_annot.rds")

# TCGA_GSVA
TCGA_gsva <- readRDS("../Zenodo_download/TCGA_GSVA_scores_safely.rds")

# Extract the GSVA scores from the safely run output
TCGA_GSVA <- TCGA_gsva[tcga_names]
TCGA_GSVA <- map(TCGA_GSVA, function(x){x$result})

# Selected enzymes of the Tryptophan degradation pathway
trp_enz_sel <- c("IL4I1","IDO1","IDO2","TDO2","TPH1","TPH2","DDC")

# the AHR signature genes
overlapping_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", stringsAsFactors = F)

#########################
## Pearson correlation ##
#########################
# Correlate between the modules and the GSVA scores (+ve & -ve correlations)
# Correlation  is  an  effect  size  and  so  we  can  verbally  describe  the  strength  of  the 
# correlation using the guide that Evans (1996) suggests for the absolute value of r:
# .00 - .19 very weak
# .20 - .39 weak
# .40 - .59 moderate
# .60 - .79 strong
# .80 - 1.0 very strong
# Evans, J. D. (1996). Straightforward statistics for the behavioral  sciences. Pacific Grove, CA: Brooks/Cole Publishing
# We will select the value of (+/-)0.2 as a cutoff with a pvalue=0.05
TCGA_MEs_GSVA_cors <- map2(TCGA_toms, TCGA_GSVA, function(tom, gsva,gene_name,cor_cutoff,cor_pvalue){
  cor_mat <- Hmisc::rcorr(as.matrix(tom$MEs), as.numeric(gsva[1,]))
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      modules = rownames(cormat)[row(cormat)[ut]],
      cor_genes = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut], stringsAsFactors = FALSE)
  }
  cor_df <- flattenCorrMatrix(cor_mat$r, cor_mat$P)
  cor_df$cor_genes[which(cor_df$cor_genes == "y")] <- gene_name
  cor_df <- cor_df[cor_df$cor_genes==gene_name,] %>%
    .[which(.$cor>=cor_cutoff | .$cor<=(cor_cutoff*-1)),] %>%
    .[.$p <=cor_pvalue,] %>% .[order(.$cor, decreasing = TRUE),]
  cor_df$modules <- gsub("^ME","",cor_df$modules)
  cor_df <- data.frame(cor_df, stringsAsFactors = F)
},gene_name="AHR_signature", cor_cutoff=0.2, cor_pvalue=0.05)

names(TCGA_MEs_GSVA_cors) <- tcga_names
saveRDS(TCGA_MEs_GSVA_cors, "../Results/RDS/TCGA_MEs_GSVA_cors.rds")

#################
## Global test ##
#################
## non directional global testing
TCGA_MEs_GSVA_GTs_nodir <- map2(1:32,tcga_names,function(i, p_n, gsva, toms){
  gsva_sc <- gsva[[i]][1,]
  test <- globaltest::gt(gsva_sc~1, ~., data=toms[[i]]$MEs, directional = FALSE)
  test_res <- globaltest::covariates(test, pdf = paste("../Results/GlobalTest/non_dir_covar_", p_n, sep = ""), zoom = TRUE)
  globaltest::subjects(test, pdf = paste("../Results/GlobalTest/non_dir_subjects_", p_n, sep = ""))
  test_res_ext <- globaltest::extract(test_res)
  test_leafs <- globaltest::leafNodes(test_res_ext)
  res <- cbind(test_leafs@result, test_leafs@extra)
}, gsva=TCGA_GSVA, toms=TCGA_toms)

names(TCGA_MEs_GSVA_GTs_nodir) <- tcga_names

saveRDS(TCGA_MEs_GSVA_GTs_nodir, "../Results/RDS/TCGA_MEs_GTs_nodir.rds")

## remove all tumors that have no hits
TCGA_MEs_GSVA_GTs_not_null_nodir <- map(TCGA_MEs_GSVA_GTs_nodir, function(x){
  if(dim(x)[[1]] == 0){
    NULL
  } else{
    x
  }
})

not_AHR_tumors_nodir <- names(TCGA_MEs_GSVA_GTs_not_null_nodir)[map(TCGA_MEs_GSVA_GTs_not_null_nodir,is.null)%>% unlist()]
not_AHR_tumors_idcs_nodir <- match(not_AHR_tumors_nodir, names(TCGA_MEs_GSVA_GTs_nodir))

## retain only the modules that have a significant assocoation with AHR activation p.value=0.05
TCGA_MEs_GSVA_GTs_not_null_sig05_nodir <- map(TCGA_MEs_GSVA_GTs_not_null_nodir, function(x){
  x[x[,1]<=0.05,]
})

saveRDS(TCGA_MEs_GSVA_GTs_not_null_sig05_nodir, "../Results/RDS/TCGA_MEs_GSVA_GTs_not_null_sig05_nodir.rds")

################################################
## Overlaping Pearson and Global test results ##
################################################
## Overlap the global test modules with the correlation modules
# Module names
TCGA_MEs_GSVA_GT_overlap_nodir <- map2(TCGA_MEs_GSVA_GTs_not_null_sig05_nodir, TCGA_MEs_GSVA_cors, function(x,y){
  GT_MEs <- gsub("ME","",rownames(x))
  Cor_MEs <- y$modules
  GT_MEs[match(Cor_MEs, GT_MEs)] %>% .[!is.na(.)]
})
saveRDS(TCGA_MEs_GSVA_GT_overlap_nodir , "../Results/RDS/TCGA_MEs_GSVA_GT_overlap_nodir.rds")

# Module Eigen genes
GT_GSVA_overlap_MEs <- map2(TCGA_MEs_GSVA_GT_overlap_nodir, TCGA_toms, function(x,y){
  y$MEs[,match(x, gsub("ME","",colnames(y$MEs)))]
})
saveRDS(GT_GSVA_overlap_MEs, "../Results/RDS/GT_GSVA_overlap_MEs.rds")


#################################
## Prepare for the circos plot ##
#################################

## Create a list where we add the pearson correlation values to the overlapping modules
TCGA_MEs_GSVA_cors_all_overlap_GT <- map2(TCGA_MEs_GSVA_cors, TCGA_MEs_GSVA_GT_overlap_nodir,function(a,b){
  a[a$modules %in% b,]
})

# Define the positive only correlating modules
TCGA_MEs_GSVA_cors_pos_only_overlap_GT <- map2(TCGA_MEs_GSVA_cors, TCGA_MEs_GSVA_GT_overlap_nodir,function(a,b){
  a[a$modules %in% b,] %>% .[.$cor>0,]
})

saveRDS(TCGA_MEs_GSVA_cors_pos_only_overlap_GT, "../Results/RDS/TCGA_MEs_GSVA_cors_pos_only_overlap_GT.rds")

## To create a circos plot for the IDO1 and TDO2, first create a list of modules with only positive
## correlations because we assume that IDO1 and TDO2 associate positively with AHR activation.
## The presence of the enzymes in these modules is defined.
modules_with_enzymes_pos <- map2(TCGA_modules_genes, TCGA_MEs_GSVA_cors_all_overlap_GT, function(df, cors, enz, sig){
  df$aaa_enzymes <- ifelse(!is.na(match(df$genes,enz)), "yes", NA)
  df$AHR_genes <- ifelse(!is.na(match(df$genes,sig$Gene)), "yes", NA)
  df <- df[which(!is.na(df$aaa_enzymes) | !is.na(df$AHR_genes)),]
  df <- df[df$module %in% cors$modules[cors$cor>=0],]
  df
}, enz=trp_enz_sel, sig=overlapping_genes)
names(modules_with_enzymes_pos) <- tcga_names

saveRDS(modules_with_enzymes_pos,"../Results/RDS/modules_with_enzymes_pos.rds")

## Matrices with the expression values of the genes in the different modules that are positively associated with AHR
TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules <- map(1:length(TCGA_modules_genes),function(i,mod_genes, gsva_corrs, voom_genes){
  res_cmg <- mod_genes[[i]][mod_genes[[i]]$module %in% gsva_corrs[[i]]$modules,]
  res_cmg_nums <- as.data.frame(table(res_cmg$module))
  res_cmg_nums <- res_cmg_nums[order(res_cmg_nums$Freq, decreasing = T),]
  res_cmg$module <- factor(res_cmg$module, levels = res_cmg_nums$Var1)
  res_cmg <- res_cmg[order(res_cmg$module),]
  res_cmg_E <- t(voom_genes[[i]]$E[,match(res_cmg$genes,colnames(voom_genes[[i]]$E))])
  res_cmg_all <- data.frame(res_cmg, res_cmg_E)
}, mod_genes=TCGA_modules_genes, gsva_corrs=TCGA_MEs_GSVA_cors_pos_only_overlap_GT, voom_genes=TCGA_voom)

names(TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules) <- tcga_names

saveRDS(TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules,"../Results/RDS/TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules.rds")

#######################################
## WGCNA Circos-plot representations ##
#######################################

save_paths_2 <- paste("../Results/TCGA_Circos_plots/Circos_plots_2/WGCNA_circos_IDO1_TDO2_", tcga_names, ".pdf", sep = "")

# Save the plots to disk
for(i in 1:32){
  WGCNA_Circ_plot_FUN(i, g_mats=TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules,
                      mod_w_enz=modules_with_enzymes_pos, ttl=gsub("TCGA_","",tcga_names),
                      save_path=save_paths_2, enz_vec=c("IDO1","TDO2"))
}

save_paths_7 <- paste("../Results/TCGA_Circos_plots/Circos_plots_7/WGCNA_circos_7_", tcga_names, ".pdf", sep = "")

# Save the plots to disk
for(i in 1:32){
  WGCNA_Circ_plot_FUN(i, g_mats=TCGA_MEs_GSVA_cors_genes_mats_pos_only_modules,
                      mod_w_enz=modules_with_enzymes_pos, ttl=gsub("TCGA_","",tcga_names),
                      save_path=save_paths_7, enz_vec=trp_enz_sel)
}
