# This script describes the generation of the immune infiltration GSVA scores for the different cell types / tumors

## Load libraries
library(GSVA)

## TCGA_voom
TCGA_voom <- readRDS("./RDS/TCGA_DGE_voom_annot.rds")

# Read the infiltration gene lists
cell_types_metagene_lists <- readRDS("./RDS/IPS_cell_types_gene_lists.rds")

## estimate the GSVA scores for the different cell populations
GSVA_voom_Immune_cells <- map(TCGA_voom, function(a,b,c){
  glists <- map(b, function(x)x$metagene)
  gsva(expr = t(a$E),glists, parallel.sz=c,method="gsva")
},b=cell_types_metagene_lists,c=6)

saveRDS(GSVA_voom_Immune_cells,"./RDS/TCGA_GSVA_voom_28_cell_types_IPS.rds")
