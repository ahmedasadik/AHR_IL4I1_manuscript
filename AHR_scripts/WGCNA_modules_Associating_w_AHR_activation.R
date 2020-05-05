#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes how modules associating with AHR activation are defined
## This is achieved by testing the association of AHR activation (GSVA scores) with the WGCNA modules.
## To overcome bias, we will try two different methods and take the overlap between them.
#################################################

## Load libraries
library(purrr)
library(globaltest)

## Source the functions and parameters files
source("./functions_and_parameters.R")

## Data
# TCGA_modules
TCGA_toms <- readRDS("./Zenodo_download/TCGA_TOMS_bicor.rds")

# TCGA_voom
TCGA_voom <- readRDS("./Zenodo_download/TCGA_DGE_voom_annots.rds")

# TCGA_GSVA
TCGA_gsva <- readRDS("./Zenodo_download/TCGA_GSVA_scores_safely.rds")
# Extract the GSVA scores from the safely run output
TCGA_GSVA <- TCGA_gsva[tcga_names]
TCGA_GSVA <- map(TCGA_GSVA, function(x){x$result})


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

saveRDS(TCGA_MEs_GSVA_cors, "./Results/RDS/TCGA_MEs_GSVA_cors.rds")

#################
## Global test ##
#################
## non directional global testing
TCGA_MEs_GSVA_GTs_nodir <- map2(1:32,names(TCGA_toms),function(i, p_n, gsva, toms){
  gsva_sc <- gsva[[i]][1,]
  test <- gt(gsva_sc~1, ~., data=toms[[i]]$MEs, directional = FALSE)
  test_res <- covariates(test, pdf = paste("./Results/GlobalTest/non_dir_covar_", p_n, sep = ""), zoom = TRUE)
  subjects(test, pdf = paste("./Results/GlobalTest/non_dir_subjects_", p_n, sep = ""))
  test_res_ext <- extract(test_res)
  test_leafs <- leafNodes(test_res_ext)
  res <- cbind(test_leafs@result, test_leafs@extra)
}, gsva=TCGA_GSVA, toms=TCGA_toms)

names(TCGA_MEs_GSVA_GTs_nodir) <- tcga_names

## remove all tumors that have no hits
TCGA_MEs_GSVA_GTs_not_null_nodir <- map(TCGA_MEs_GSVA_GTs_nodir, function(x){ifelse(dim(x)[[1]] == 0,NULL,x)})

not_AHR_tumors_nodir <- names(TCGA_MEs_GSVA_GTs_not_null_nodir)[map(TCGA_MEs_GSVA_GTs_not_null_nodir,is.null)%>% unlist()]
not_AHR_tumors_idcs_nodir <- match(not_AHR_tumors_nodir, names(TCGA_MEs_GSVA_GTs_nodir))

## retain only the modules that have a significant assocoation with AHR activation p.value=0.05
TCGA_MEs_GSVA_GTs_not_null_sig05_nodir <- map(TCGA_MEs_GSVA_GTs_not_null_nodir, function(x){
  x[x[,1]<=0.05,]
})

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
saveRDS(TCGA_MEs_GSVA_GT_overlap_nodir , "./Results/RDS/TCGA_MEs_GSVA_GT_overlap_nodir.rds")

# Module Eigen genes
GT_GSVA_overlap_MEs <- map2(TCGA_MEs_GSVA_GT_overlap_nodir, TCGA_toms, function(x,y){
  y$MEs[,match(x, gsub("ME","",colnames(y$MEs)))]
})
saveRDS(GT_GSVA_overlap_MEs, "./Results/RDS/GT_GSVA_overlap_MEs.rds")

