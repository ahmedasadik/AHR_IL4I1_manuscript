# This script describes the genration of the ontology terms of the AHR associating modules

## Load libraries
library(purrr)
library(parallel)
library(org.Hs.eg.db)
library(clusterProfiler)

# Annotation update function
annot_query_FUN <- function(gS, sym, eid, syn, psym=NA){
  if (length(grep(paste("^(", gS, ")$", sep = ""), sym)) > 0){
    idx <- grep(paste("^(", gS, ")$", sep = ""), sym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("^(", gS, ")$", sep = ""), psym)) > 0){
    idx <- grep(paste("^(", gS, ")$", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("^(", gS, ")\\|", sep = ""), psym)) > 0){
    idx <- grep(paste("\\|(", gS, ")\\|", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("\\|(", gS, ")\\|", sep = ""), psym)) > 0){
    idx <- grep(paste("\\|(", gS, ")\\|", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (!is.na(psym) && length(grep(paste("\\|(", gS, ")$", sep = ""), psym)) > 0){
    idx <- grep(paste("\\|(", gS, ")$", sep = ""), psym)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else if (length(grep(paste("^(", gS, ")$", sep = ""), syn)) > 0){
    idx <- grep(paste("^(", gS, ")$", sep = ""), syn)
    return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
  } else {
    idx <- grep(paste("^(", gS, ")\\|", sep = ""), syn)
    idx2 <- grep(paste("\\|(", gS, ")\\|", sep = ""), syn)
    idx3 <- grep(paste("\\|(", gS, ")$", sep = ""), syn)
    if (length(idx)>0){
      return(c(paste(sym[idx], collapse = ";"), paste(eid[idx], collapse = ";")))
    } else if (length(idx2)>0) {
      return(c(paste(sym[idx2], collapse = ";"), paste(eid[idx2], collapse = ";")))
    } else {
      return(c(paste(sym[idx3], collapse = ";"), paste(eid[idx3], collapse = ";")))
    }
  }
}

# Module Eigen genes
GT_GSVA_overlap_MEs <- readRDS("./RDS/GT_GSVA_overlap_MEs.rds")

# Gene names in the different modules
TCGA_modules_genes <- readRDS("./RDS/TCGA_modules_genes.rds")

# HGNC
hgnc <- read.delim("/home/data/Resources/HGNC/hgnc_complete_set_04122017.txt")

#########################################################################
## Generate the GO ontologies that are over-represented in each module ##
#########################################################################

## Generate a list for each tumor with the genes in each module and add their EIDs
modules_genes_lists_separate <- map2(GT_GSVA_overlap_MEs, TCGA_modules_genes, function(a1,a2){
  MEs <- names(a1) %>% gsub("ME","",.)
  ME_gene_lists <- purrr::map(MEs,function(x){a2[a2$module==x,]})
  names(ME_gene_lists) <- MEs
  ME_gene_lists
})
Tumor_names_reps_GO_lists <- rep(names(modules_genes_lists_separate), unlist(lapply(modules_genes_lists_separate, length)))
Tumor_module_names_GO_lists <- unlist(lapply(modules_genes_lists_separate, names))
module_genes_for_GO <- map2(Tumor_names_reps_GO_lists, Tumor_module_names_GO_lists, function(a1, a2){
  modules_genes_lists_separate[[a1]][[a2]]$genes})
names(module_genes_for_GO) <- paste(Tumor_names_reps_GO_lists, Tumor_module_names_GO_lists, sep = "_")

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl,varlist = c("hgnc","annot_query_FUN"))
EID_annot <- parLapply(cl, module_genes_for_GO, fun= function(gsl, sym, eid, syn, psym){
  t(sapply(gsl,annot_query_FUN, sym = sym, eid = eid, syn = syn, psym = psym))},
  sym = hgnc$symbol, eid = hgnc$entrez_id, syn = hgnc$alias_symbol,
  psym = hgnc$prev_symbol)
stopCluster(cl)

EID_annot_na <- purrr::map(EID_annot,function(a1){
  empty_as_na <- function(x){
    if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
    ifelse(as.character(x)!="", x, NA)
  }
  df <- as.data.frame(a1) %>% mutate_all(funs(empty_as_na)) %>% .[!is.na(.[,2]),]
  colnames(df) <- c("hGene", "EID")
  df
})
saveRDS(EID_annot_na, "./RDS/EID_annot_na.rds")

## Run gene ontology analysis on the EIDs for the different modules
GO_ALL_tumor_modules <- lapply(EID_annot_na, safely(function(a1){
  GO_res <- enrichGO(gene = a1$EID, OrgDb=org.Hs.eg.db, ont = "BP", pAdjustMethod = "bonferroni", readable = TRUE, pvalueCutoff = 0.01)
  GO_res_simple <- clusterProfiler::simplify(GO_res)
  GO_res_filtered <- gofilter(GO_res_simple)
  print("Done")
  res <- list(BP=GO_res, BP_simple=GO_res_simple, BP_sim_fil=GO_res_filtered)
}))
saveRDS(GO_ALL_tumor_modules, "./RDS/GO_all_tumor_modules.rds")

# Saving all hits in separate files before any filters are applied
GO_BP_dfs <- map(GO_ALL_tumor_modules, safely(function(x)x$result$BP@result))
GO_BP_dfs_idx <- map(GO_BP_dfs, function(x)is.null(x$error))%>% unlist()
GO_BP_dfs <- GO_BP_dfs[GO_BP_dfs_idx]

walk2(GO_BP_dfs,paste("./BP/",names(GO_BP_dfs),".csv",sep = ""), function(a1,a2){
  write.csv(as.data.frame(a1$result),a2)})

# Saving all hits in separate files after merging and filtering are applied
GO_BP_sim_dfs <- map(GO_ALL_tumor_modules, safely(function(x)x$result$BP_sim_fil@result))
GO_BP_sim_dfs_idx <- map(GO_BP_sim_dfs, function(x)is.null(x$error))%>% unlist()
GO_BP_sim_dfs <- GO_BP_sim_dfs[GO_BP_sim_dfs_idx]
walk2(GO_BP_sim_dfs,paste("./BP_sim/",names(GO_BP_sim_dfs),".csv",sep = ""), function(a1,a2){
  write.csv(as.data.frame(a1$result),a2)})
