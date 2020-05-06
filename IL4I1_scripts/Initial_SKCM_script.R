# This script describes a full run for metastatic melanoma similar to that of the the TP TCGA cohort

## Load libraries
library(purrr)
library(edgeR)
library(WGCNA)
library(ggpubr)

## colors ####
cluster_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                    "#e377c2","#7f7f7f","#bcbd22","#17becf")

## load TCGA_voom
load("/home/analyses/Projects/PRMT5/tcga_results/tcga_rda/TCGA-SKCMTranscriptome_ProfilingWed_Oct__4_161844_2017.RData")

TCGA_SKCM_colData <- SummarizedExperiment::colData(data)
TCGA_SKCM_genes <- SummarizedExperiment::rowData(data)
TCGA_SKCM_counts <- SummarizedExperiment::assay(data)

## Extracting only metastatic patients ####
samplesTM <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(TCGA_SKCM_counts),typesample = c("TM")) 
counts_TM <- TCGA_SKCM_counts[,match(samplesTM, colnames(TCGA_SKCM_counts))]
TCGA_SKCM_colData_TM <- TCGA_SKCM_colData[match(samplesTM, TCGA_SKCM_colData$barcode),]

skcm_dge <- DGEList(counts_TM, samples = TCGA_SKCM_colData_TM, genes = TCGA_SKCM_genes)
saveRDS(skcm_dge,"../SKCM/SKCM_DGE.rds")

## Generate the AHR GSVA scores
AHR_genes <- read.delim("/home/analyses/Projects_analyses/AHR_IL4I1/AHR/Publication_results/Resources/Signature/overlapping_AHR_signature.txt", stringsAsFactors = F)
AHR_ensgs <- skcm_dge$genes[skcm_dge$genes$external_gene_name %in% AHR_genes$Gene,]

library(GSVA)
AHR_gsva_met <- gsva(expr = cpm(skcm_dge$counts, prior.count = 1, log = T), list(AHR_ensgs$ensembl_gene_id), method="gsva", parallel.sz=6)
saveRDS(AHR_gsva_met,"../SKCM/SKCM_GSVA.rds")

## Correlating the AHR score with selected trp enzymes ##
trp_enz_sel <- c("IL4I1","IDO1","TDO2","KMO","KYNU","MAOB","MAOA","CAT","TPH1","TPH2","DDC","KYAT3","KYAT1","AADAT","AFMID")
Enz_ensg_vec <- TCGA_SKCM_genes[TCGA_SKCM_genes$external_gene_name %in% trp_enz_sel,]

skcm_dge2 <-skcm_dge 
skcm_dge2$counts <- cpm(skcm_dge2$counts, prior.count = 1, log = T)
Enz_counts_matrix <- skcm_dge2$counts[rownames(skcm_dge2$counts) %in% Enz_ensg_vec$ensembl_gene_id,]
rownames(Enz_counts_matrix) <- Enz_ensg_vec$external_gene_name

df <- data.frame(AHR_score=AHR_gsva_met[1,], t(Enz_counts_matrix), stringsAsFactors = F)
rownames(df) <- substring(rownames(df), 1,12)

df_corrs <- Hmisc::rcorr(as.matrix(df[,-1]), as.numeric(df[,1]))

df_corrs_df <- data.frame(Enzymes=rownames(df_corrs$r)[1:16],corrs=df_corrs$r[,16], pvalues=df_corrs$P[,16], stringsAsFactors = F)

df_corrs_df <- df_corrs_df[-16,]
df_corrs_df <- df_corrs_df[order(df_corrs_df$corrs, decreasing = T),]
df_corrs_df$corrs <- ifelse(df_corrs_df$pvalues<0.05, df_corrs_df$corrs, 0)
df_corrs_df$fill_colors <- ifelse(df_corrs_df$corrs > 0, "1","2")

ggbarplot(df_corrs_df, "Enzymes", "corrs", fill = "fill_colors", palette = "npg", ylab = "Correlation with AHR Activation\n",
          xlab = FALSE, legend="none", x.text.angle = 90)

## WGCNA prep ####
keep <- filterByExpr(skcm_dge)
dge <- skcm_dge[keep,,keep.lib.size=FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, plot=FALSE)

rownames(v$E) <- v$genes$external_gene_name [match(rownames(v$E), v$genes$ensembl_gene_id)]
v$E <- t(v$E)

saveRDS(v, "../SKCM/SKCM_voom.rds")

# Enabling WGCNA threading
enableWGCNAThreads(nThreads = 4)
# Estimating soft threshold power
# Selecting the powers for soft thresholding
powers <- c(c(1:10), seq(from = 12, to=40, by=2))

# Select soft threshold
TCGA_sft <- pickSoftThreshold(v$E, dataIsExpr = TRUE, networkType = "signed hybrid",
                    blockSize = dim(v$E)[[2]], powerVector = powers, verbose = 5)

saveRDS(TCGA_sft, "../SKCM/SKCM_sft.rds")

## RUN WGCNA and generate coexpression modules ####
TCGA_TOMS <- blockwiseModules(datExpr = (v$E), corType = "bicor",
                   maxBlockSize = dim(v$E)[[2]], nThreads = 4,
                   randomSeed = 0861, power = TCGA_sft$powerEstimate,
                   networkType = "signed hybrid",
                   TOMType = "signed", TOMDenom = "mean", saveTOMs = TRUE,
                   saveTOMFileBase = "sckm_metastasis",
                   pamStage = TRUE, pamRespectsDendro = FALSE,verbose = 3)

saveRDS(TCGA_TOMS,"../SKCM/SKCM_TOMS_bicor.rds")

## Extract the gene names for the different modules ####
SKCM_modules_genes <- data.frame(genes=colnames(v$E), module=TCGA_TOMS$colors, stringsAsFactors = F)
saveRDS(SKCM_modules_genes, "../SKCM/SKCM_modules_genes.rds")

## Pearson correlations of the modules and AHR
SKCM_MEs_GSVA_cors <- map2(list(TCGA_TOMS), list(AHR_gsva_met), function(tom, gsva,gene_name,cor_cutoff,cor_pvalue)
{
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
},gene_name="AHR_signature",cor_cutoff=0.2,cor_pvalue=0.05)
names(SKCM_MEs_GSVA_cors) <- "SKCM"

saveRDS(SKCM_MEs_GSVA_cors, "../SKCM/SKCM_MEs_GSVA_cors.rds")

library(globaltest)

## non directional
SKCM_MEs_GSVA_GTs_nodir <- map2(list(1),list("SKCM"),function(i, p_n, gsva, toms){
  gsva_sc <- gsva[1,]
  test <- gt(gsva_sc~1, ~., data=toms$MEs, directional = FALSE)
  test_res <- covariates(test, pdf = paste("../SKCM/non_dir_covar_", p_n, sep = ""), zoom = TRUE)
  subjects(test, pdf = paste("../SKCM/non_dir_subjects_", p_n, sep = ""))
  test_res_ext <- extract(test_res)
  test_leafs <- leafNodes(test_res_ext)
  res <- cbind(test_leafs@result, test_leafs@extra)
}, gsva=AHR_gsva_met, toms=TCGA_TOMS)
names(SKCM_MEs_GSVA_GTs_nodir) <- "SKCM"

## remove all tumors that have no hits
SKCM_MEs_GSVA_GTs_not_null_nodir <- map(SKCM_MEs_GSVA_GTs_nodir, function(x){
  if(dim(x)[[1]] == 0){
    NULL
  } else{
    x
  }
})

not_AHR_tumors_nodir <- names(SKCM_MEs_GSVA_GTs_not_null_nodir)[map(SKCM_MEs_GSVA_GTs_not_null_nodir,is.null)%>% unlist()]
not_AHR_tumors_idcs_nodir <- match(not_AHR_tumors_nodir, names(SKCM_MEs_GSVA_GTs_nodir))

## retain only the modules that have a significant assocoation with AHR activation p.value=0.05
SKCM_MEs_GSVA_GTs_not_null_sig05_nodir <- map(SKCM_MEs_GSVA_GTs_not_null_nodir, function(x){
  x[x[,1]<=0.05,]
})

## Overlap the global test modules with the correlation modules ####
SKCM_MEs_GSVA_GT_overlap_nodir <- map2(SKCM_MEs_GSVA_GTs_not_null_sig05_nodir, SKCM_MEs_GSVA_cors, function(x,y){
  GT_MEs <- gsub("ME","",rownames(x))
  Cor_MEs <- y$modules
  GT_MEs[match(Cor_MEs, GT_MEs)] %>% .[!is.na(.)]
})

GT_GSVA_overlap_MEs <- map2(SKCM_MEs_GSVA_GT_overlap_nodir,list(TCGA_TOMS), function(x,y){
  y$MEs[,match(x, gsub("ME","",colnames(y$MEs)))]
})

saveRDS(SKCM_MEs_GSVA_GT_overlap_nodir , "../SKCM/SKCM_MEs_GSVA_GT_overlap_nodir.rds")
saveRDS(GT_GSVA_overlap_MEs, "../SKCM/SKCM_GT_GSVA_overlap_MEs.rds")

## Create a list of overlapping modules and their correlations
SKCM_MEs_GSVA_cors_all_overlap_GT <- map2(SKCM_MEs_GSVA_cors,SKCM_MEs_GSVA_GT_overlap_nodir,function(a,b){
  a[a$modules %in% b,]
})

SKCM_MEs_GSVA_cors_pos_only_overlap_GT <- map2(SKCM_MEs_GSVA_cors,SKCM_MEs_GSVA_GT_overlap_nodir,function(a,b){
  a[a$modules %in% b,] %>% .[.$cor>0,]
})

modules_with_enzymes_pos <- map2(list(SKCM_modules_genes), SKCM_MEs_GSVA_cors_all_overlap_GT, function(df, cors, enz, sig){
  df$aaa_enzymes <- ifelse(!is.na(match(df$genes,enz)), "yes", NA)
  df$AHR_genes <- ifelse(!is.na(match(df$genes,sig$Gene)), "yes", NA)
  df <- df[which(!is.na(df$aaa_enzymes) | !is.na(df$AHR_genes)),]
  df <- df[df$module %in% cors$modules[cors$cor>=0],]
  df
}, enz=trp_enz_sel, sig=AHR_genes)

#### WGCNA Circos-plot representations ####
save_paths <- paste("../SKCM/WGCNA_circos_", "SKCM", ".tiff", sep = "")

## Matrices with the expression values of the genes in the different modules that are positively associated with AHR
SKCM_MEs_GSVA_cors_genes_mats_pos_only_modules <- map(1,function(i,mod_genes, gsva_corrs, voom_genes){
  res_cmg <- mod_genes[mod_genes$module %in% gsva_corrs[[i]]$modules,]
  res_cmg_nums <- as.data.frame(table(res_cmg$module))
  res_cmg_nums <- res_cmg_nums[order(res_cmg_nums$Freq, decreasing = T),]
  res_cmg$module <- factor(res_cmg$module, levels = res_cmg_nums$Var1)
  res_cmg <- res_cmg[order(res_cmg$module),]
  res_cmg_E <- t(voom_genes[[i]]$E[,match(res_cmg$genes,colnames(voom_genes[[i]]$E))])
  res_cmg_all <- data.frame(res_cmg, res_cmg_E)
},
mod_genes=SKCM_modules_genes, gsva_corrs=SKCM_MEs_GSVA_cors_pos_only_overlap_GT, voom_genes=list(v))
names(SKCM_MEs_GSVA_cors_genes_mats_pos_only_modules) <- "SKCM"

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

library(circlize)
library(ComplexHeatmap)
WGCNA_Circ_plot_FUN(1, g_mats=SKCM_MEs_GSVA_cors_genes_mats_pos_only_modules,
                    mod_w_enz=modules_with_enzymes_pos, ttl="SKCM", save_path=save_paths, enz_vec=trp_enz_sel)

## Consensus clustering
library(ConsensusClusterPlus)
library(parallel)

ccp_fun_km_dist <- function(mods, dirs){
  d_m <- 1-cor(as.matrix(mods)) %>% as.dist()
  ccp_run <- ConsensusClusterPlus(d_m, maxK = 10, reps = 1000,
                                  clusterAlg = "kmdist", writeTable = FALSE, seed = 0860,
                                  title = paste(dirs,"_ccp", sep = ""), plot="pdf")
  calcIC_run <- calcICL(ccp_run,plot = "pdf", title = paste(dirs,"_calcIC", sep = ""))
  res <- list(ccp=ccp_run, calcIC=calcIC_run)
}

ccp_pos_dist <- ccp_fun_km_dist(t(GT_GSVA_overlap_MEs$SKCM),"../SKCM/SKCM")

saveRDS(ccp_pos_dist, "../SKCM/SKCM_ccp_pos_dist.rds")

## SKCM is 4 groups
SKCM_clusters <- paste("K",c(4), sep = "")

## PCA validation
library(FactoMineR)
library(factoextra)
library(RColorBrewer)

set.seed(0860)
PCA_res <- purrr::map(GT_GSVA_overlap_MEs, PCA)

### PCA scree plots ####
PCA_Scree_plots <- purrr::map2(PCA_res,names(PCA_res), function(a,b){
  fviz_screeplot(a,main=b, addlabels = TRUE, linecolor="red", barcolor="grey")
})
purrr::walk2(paste("../SKCM/Plots_scree_",names(PCA_Scree_plots),".png", sep = ""), PCA_Scree_plots,
             function(a,b){
               png(a, res = 600, width = 12, height = 8, units = "in")
               print(b)
               dev.off()
             })
### PCA paired plots ####
PCA_pairs_plots <- purrr::map(1,function(idcs,a1, a2, a3){
  df <- as.data.frame(a1[[idcs]]$ind$coord)
  ccp_input <- a2[[idcs]]$ccp
  ccp_idcs <- cbind(2,3,4,5)
  ccp_idcs_mat <- apply(ccp_idcs,2, function(x){
    ccp_input[[x]]$consensusClass
  })
  paired_plots <- apply(ccp_idcs_mat,2, function(ccp_i){
    p <- GGally::ggpairs(df,aes(col=as.factor(ccp_i)),title =a3[idcs])
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
},a1=PCA_res,a2=list(ccp_pos_dist),a3="SKCM")
names(PCA_pairs_plots) <- names(PCA_res)

TCGA_PCA_tumor_idcs_for_saving <- rep(c(1),each=4)
PCA_pairs_plots_idcs_for_saving <- c(1:4)
pairs_plots_to_save <- purrr::map2(TCGA_PCA_tumor_idcs_for_saving,PCA_pairs_plots_idcs_for_saving,
                                   function(a,b){
                                     p <- PCA_pairs_plots[[a]][b]
                                   })

to_save_names <- purrr::map2(TCGA_PCA_tumor_idcs_for_saving,
                             PCA_pairs_plots_idcs_for_saving,
                             function(a,b){
                               paste(names(PCA_res)[a], names(PCA_pairs_plots[[a]])[b], sep = "_")
                             }) %>% unlist()

purrr::walk2(paste("../SKCM/Plots_pairs_",to_save_names,".png", sep = ""),
             pairs_plots_to_save,
             function(a,b){
               png(a, res = 600, width = 12, height = 8, units = "in")
               print(b)
               dev.off()
             })

### PCA 3D plots ####
require(plotly)
PCA_3D_plots <- purrr::map(1,function(idcs,a1, a2, a3){
  df <- as.data.frame(a1[[idcs]]$ind$coord)
  ccp_input <- a2[[idcs]]$ccp
  ccp_idcs <- cbind(2,3,4,5)
  ccp_idcs_mat <- apply(ccp_idcs,2, function(x){
    ccp_input[[x]]$consensusClass
  })
  paired_plots <- apply(ccp_idcs_mat,2, function(ccp_i){
    p <- plot_ly(df, split=ccp_i, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, type = "scatter3d", mode = "lines")%>% layout(title = a3[[idcs]]) 
  }) 
  names(paired_plots) <- c("K2", "K3", "K4", "K5")
  paired_plots
},a1=PCA_res,a2=list(ccp_pos_dist), a3="SKCM")
names(PCA_3D_plots) <- "SKCM"

### The number of clusters of each TCGA tumor ####
### Run DGE ####
RNA_DGE_WGCNA_voomed_data_FUN <- function(des_mat, dge){
  if(length(colnames(des_mat))==2){
    cont_mat <- makeContrasts(cls_1-cls_2, levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit, robust = TRUE)
    tt <- topTable(efit, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt, eBFit=efit)
  } else if(length(colnames(des_mat))==3){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_2-cls_3, levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit)
    tt1 <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "M")
    tt2 <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")
    tt3 <- topTable(efit, coef = 3, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt1, cnt2=tt2, cnt3=tt3, eBFit=efit)
  } else if(length(colnames(des_mat))==4){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_1-cls_4, cls_2-cls_3, cls_2-cls_4, cls_3-cls_4,levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit)
    tt1 <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "M")
    tt2 <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")
    tt3 <- topTable(efit, coef = 3, number = Inf, adjust.method = "BH", sort.by = "M")
    tt4 <- topTable(efit, coef = 4, number = Inf, adjust.method = "BH", sort.by = "M")
    tt5 <- topTable(efit, coef = 5, number = Inf, adjust.method = "BH", sort.by = "M")
    tt6 <- topTable(efit, coef = 6, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt1, cnt2=tt2, cnt3=tt3, cnt4=tt4, cnt5=tt5, cnt6=tt6, eBFit=efit)
  } else if(length(colnames(des_mat))==5){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_1-cls_4,cls_1-cls_5,
                              cls_2-cls_3, cls_2-cls_4,cls_2-cls_5, cls_3-cls_4,
                              cls_3-cls_5,cls_4-cls_5,levels = des_mat)
    vfit<-lmFit(t(dge$E), des_mat)
    vfit<-contrasts.fit(vfit,contrasts=cont_mat)
    efit<-eBayes(vfit)
    tt1 <- topTable(efit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "M")
    tt2 <- topTable(efit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "M")
    tt3 <- topTable(efit, coef = 3, number = Inf, adjust.method = "BH", sort.by = "M")
    tt4 <- topTable(efit, coef = 4, number = Inf, adjust.method = "BH", sort.by = "M")
    tt5 <- topTable(efit, coef = 5, number = Inf, adjust.method = "BH", sort.by = "M")
    tt6 <- topTable(efit, coef = 6, number = Inf, adjust.method = "BH", sort.by = "M")
    tt7 <- topTable(efit, coef = 7, number = Inf, adjust.method = "BH", sort.by = "M")
    tt8 <- topTable(efit, coef = 8, number = Inf, adjust.method = "BH", sort.by = "M")
    tt9 <- topTable(efit, coef = 9, number = Inf, adjust.method = "BH", sort.by = "M")
    tt10 <- topTable(efit, coef =10, number = Inf, adjust.method = "BH", sort.by = "M")
    tts <- list(cnt1=tt1, cnt2=tt2, cnt3=tt3, cnt4=tt4, cnt5=tt5, cnt6=tt6,
                cnt7=tt7, cnt8=tt8, cnt9=tt9, cnt10=tt10,
                eBFit=efit)
  }
}
SKCM_DGE_ccp <- purrr::map2(list(v), list(ccp_pos_dist), function(a1,a2){
  ccp_input <- a2$ccp
  ccp_idcs <- cbind(2,3,4,5)
  ccp_idcs_mat <- apply(ccp_idcs,2, function(x){
    ccp_input[[x]]$consensusClass
  })
  des_mats <- apply(ccp_idcs_mat,2,function(y){
    des<-model.matrix(~0+as.factor(y))
    colnames(des) <- paste("cls",levels(as.factor(y)),sep="_")
    des})
  names(des_mats) <- c("K2", "K3", "K4", "K5")
  tts <- lapply(des_mats, RNA_DGE_WGCNA_voomed_data_FUN, dge=a1)
})
names(SKCM_DGE_ccp) <- "SKCM"
saveRDS(TCGA_DGE_ccp, "./RDS/TCGA_DGE_ccp.rds")

### selected TCGA_clusters comparisons ####
SKCM_DGE_ccp_clusters <- map2(SKCM_DGE_ccp,list(SKCM_clusters), function(a1,a2){
  tts <- a1[[a2]] %>% .[-length(.)]
})

saveRDS(SKCM_DGE_ccp_clusters, "../SKCM/SKCM_DGE_ccp_clusters.rds")

walk2(SKCM_DGE_ccp_clusters, names(SKCM_DGE_ccp_clusters), function(a1,a2){
  walk2(a1, paste("../SKCM/",a2,"_tt_", names(a1),".csv",sep = ""), write.csv)
})

## Boxplot the GSVA scores distributed as per the new selected cluster assignments ####
ccp_cluster_assignments <- map2(list(ccp_pos_dist),list(4), function(a1, a2){a1$ccp[[a2]]$consensusClass})

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

GSVA_clusters_plots <- map(1, function(i,a1,a2,a3){
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
},
a1=list(AHR_gsva_met), a2=ccp_cluster_assignments, a3="SKCM")
names(GSVA_clusters_plots) <- "SKCM"

walk2(GSVA_clusters_plots, paste("../SKCM/SKCM_GSVA_boxplots_ccp_clust_K4", ".png", sep = ""),
      function(a1,a2){
        png(a2, res = 600, width = 12, height = 8, units = "in")
        print(a1)
        dev.off()
      })

## Prepare for running roast ####
## Create design matrices
TCGA_des_mats <- map2(list(ccp_pos_dist), list(SKCM_clusters), function(a1,a2){
  idx <- gsub("K","",a2) %>% as.numeric()
  clss <- a1$ccp[[idx]]$consensusClass
  des<-model.matrix(~0+as.factor(clss))
  colnames(des) <- paste("cls",levels(as.factor(clss)),sep="_")
  des})

## Create contrast matrices
TCGA_cont_mats <- map(TCGA_des_mats, function(des_mat){
  if(length(colnames(des_mat))==2){
    cont_mat <- makeContrasts(cls_1-cls_2, levels = des_mat)
  } else if(length(colnames(des_mat))==3){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_2-cls_3, levels = des_mat)
  } else if(length(colnames(des_mat))==4){
    cont_mat <- makeContrasts(cls_1-cls_2, cls_1-cls_3, cls_1-cls_4,
                              cls_2-cls_3, cls_2-cls_4, cls_3-cls_4,levels = des_mat)
  }
})

## Run roast for the selected clusters ####
TCGA_roast_FUN <- function(idx, v_eset, glist, des_mat, cont_mat){
  index.vector <- colnames(v_eset$E) %in% glist
  AHR_roast <- roast(t(v_eset$E), index = index.vector, design = des_mat,
                     contrast = cont_mat[,idx], set.statistic = "floormean")
  roast_cols <- c("NGenes",	"PropDown",	"PropUp",	"Direction",	"PValue",	"FDR",	"PValue.Mixed",	"FDR.Mixed")
  AHR_roast_df <- as.matrix(rep(0,length(roast_cols)),nrow=1) %>% t() %>% as.data.frame()
  colnames(AHR_roast_df) <- roast_cols
  AHR_roast_df$NGenes <- AHR_roast$ngenes.in.set
  AHR_roast_df$PropDown <- AHR_roast$p.value$Active.Prop[1] %>% round(.,4)
  AHR_roast_df$PropUp <- AHR_roast$p.value$Active.Prop[2] %>% round(.,4)
  AHR_roast_df$Direction <- ifelse(AHR_roast_df$PropDown > AHR_roast_df$PropUp, "Down", "Up")
  AHR_roast_df$PValue <- AHR_roast_df$FDR <- ifelse(AHR_roast_df$Direction=="Down",
                                                    round(AHR_roast$p.value$P.Value[3], 4),
                                                    round(AHR_roast$p.value$P.Value[3],4))
  AHR_roast_df$PValue.Mixed <- AHR_roast_df$FDR.Mixed <- round(AHR_roast$p.value$P.Value[4],4)
  AHR_roast_df 
}

TCGA_roast_res <- map(1, function(i, v_eset, glist, des_mat, cont_mat){
  if(length(colnames(cont_mat[[i]])) == 1){
    cnts <- TCGA_roast_FUN(idx = 1, v_eset = v_eset[[i]], glist = glist, des_mat = des_mat[[i]], cont_mat[[i]])  
  } else {
    cnts <- map(1:length(colnames(cont_mat[[i]])),TCGA_roast_FUN,
                v_eset = v_eset[[i]], glist = glist, des_mat = des_mat[[i]], cont_mat[[i]]) 
  }
  cnts
},
v_eset=list(v),glist=AHR_genes$Gene, des_mat=TCGA_des_mats, cont_mat=TCGA_cont_mats)

names(TCGA_roast_res) <- "SKCM"

cnts_names <- map2(SKCM_DGE_ccp_clusters, names(SKCM_DGE_ccp_clusters), function(a1,a2){
  cnts_names <- names(a1)
  paste(a2, cnts_names, sep = "_")
}) %>% unlist()

library(dplyr)
TCGA_roast_res_df <- do.call(bind_rows, TCGA_roast_res) %>% data.frame(Contrasts=cnts_names,.,stringsAsFactors=FALSE)
TCGA_roast_res_df$pct_diff <- (abs(TCGA_roast_res_df$PropDown-TCGA_roast_res_df$PropUp)/(TCGA_roast_res_df$PropDown+TCGA_roast_res_df$PropUp))*100

write.csv(TCGA_roast_res_df, "../SKCM/SKCM_roast_results.csv")

## TCGA_barcode_plots for the selected clusters that were roasted ####
walk(1, function(i, tt_list, tcga_names, glist){
  if(length(names(tt_list[[i]])) == 1){
    index.vector2 <- tt_list[[i]][[1]]$ID %in% glist
    png(paste("../SKCM/barcodeplots_", tcga_names[i],"_cnt1.png", sep=""),width = 12, height = 8, units = "in", res = 600)
    barcodeplot(tt_list[[i]][[1]]$t, index = index.vector2,main=paste(tcga_names[i],"cnt1",sep = "_"),
                gene.weights = tt_list[[i]][[1]][index.vector2,"logFC"])
    dev.off()
  } else {
    walk(1:length(names(tt_list[[i]])),function(idx){
      index.vector2 <- tt_list[[i]][[idx]]$ID %in% glist
      png(paste("../SKCM/barcodeplots_", tcga_names[i],"_",names(tt_list[[i]][idx]),".png", sep=""),width = 12, height = 8, units = "in", res = 600)
      barcodeplot(tt_list[[i]][[idx]]$t, index = index.vector2, main=paste(tcga_names[i],names(tt_list[[i]][idx]),sep = "_"),
                  gene.weights = tt_list[[i]][[idx]][index.vector2,"logFC"])
      dev.off()
    })
  }
},
tt_list=SKCM_DGE_ccp_clusters,tcga_names="SKCM", glist=AHR_genes$Gene)
## After roast, order the groups based on AHR activation for the different tumor types ####
TCGA_DGE_ccp_clusters_AHR_activation_order <- map("SKCM", function(a1){
  df <- rep(NA,4) %>% as.data.frame() %>% t() %>% as.data.frame()
  colnames(df) <- paste("C",1:4, sep = "")
  rownames(df) <- a1
  idcs <- grep(a1, TCGA_roast_res_df$Contrasts)
  if(length(idcs)==1){
    if(TCGA_roast_res_df$Direction[idcs]=="Up"){
      df <- c(1,2,NA,NA)
    } else {
      df <- c(2,1,NA,NA)
    }
  } else if(length(idcs)==3){
    if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up"))){
      df <- c(1,2,3,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Down"))){
      df <- c(1,3,2,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Up"))){
      df <- c(2,1,3,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Up"))){
      df <- c(2,3,1,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Down"))){
      df <- c(3,1,2,NA)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down"))){
      df <- c(3,2,1,NA)
    }
  } else if(length(idcs)==6){
    if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Up","Up","Up"))){
      df <- c(1,2,3,4)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Up","Up","Down"))){
      df <- c(1,2,4,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Down","Up","Up"))){
      df <- c(1,3,2,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Down","Down","Up"))){
      df <- c(1,3,4,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Up","Down","Down"))){
      df <- c(1,4,2,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Up","Down","Down","Down"))){
      df <- c(1,4,3,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Up","Up","Up","Up"))){
      df <- c(2,1,3,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Up","Up","Up","Down"))){
      df <- c(2,1,4,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Up","Up","Up","Up"))){
      df <- c(2,3,1,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Up","Up","Up"))){
      df <- c(2,3,4,1)
    } else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Down","Up","Up","Down"))){
      df <- c(2,4,1,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Up","Up","Down"))){
      df <- c(2,4,3,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Up","Down","Up","Up"))){
      df <- c(3,1,2,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Up","Down","Down","Up"))){
      df <- c(3,1,4,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Up","Down","Up","Up"))){
      df <- c(3,2,1,4)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Down","Up","Up"))){
      df <- c(3,2,4,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Down","Down","Down","Up"))){
      df <- c(3,4,1,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Down","Down","Up"))){
      df <- c(3,4,2,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Down","Up","Down","Down"))){
      df <- c(4,1,2,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Up","Down","Down","Down","Down"))){
      df <- c(4,1,3,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Up","Down","Up","Down","Down"))){
      df <- c(4,2,1,3)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Up","Down","Down"))){
      df <- c(4,2,3,1)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Up", "Down","Down","Down","Down","Down"))){
      df <- c(4,3,1,2)
    }else if(identical(TCGA_roast_res_df$Direction[idcs],c("Down", "Down","Down","Down","Down","Down"))){
      df <- c(4,3,2,1)
    }
    df
  }
}) %>% do.call(rbind,.)
TCGA_DGE_ccp_clusters_AHR_activation_order_df <- data.frame(Tumor="SKCM", TCGA_DGE_ccp_clusters_AHR_activation_order)
write.table(TCGA_DGE_ccp_clusters_AHR_activation_order_df, "../SKCM/SKCM_DGE_ccp_clusters_AHR_activation_order_df.txt")

## Heatmaps
hm_RPPA_FUN <- function(tt_df, clust_grp, ttl, col_p, sep_vec,...){
  col_clust <- clust_grp
  ordered_patients <- tt_df[order(col_clust),]
  to_plot <- t(scale(as.matrix(ordered_patients)))
  
  gplots::heatmap.2(to_plot,col=col_p, trace = "none",Colv = NULL,main = ttl, colsep = sep_vec,
            scale = "row", cexRow =0.8, dendrogram = "none", key.title = NA,margins = c(8, 8),
            density.info = c("density"), labRow = rownames(to_plot), labCol = FALSE, 
            ColSideColors = cluster_colors[col_clust[order(col_clust)]])
}
walk(1,function(i,a1,a2,a3,a4){
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
  png(paste("../SKCM/SKCM_MEs",a3[i],a4[i],".png",sep = "_"),width = 12, height = 8, units = "in", res = 600)
  hm_RPPA_FUN(as.matrix((a1[[i]])), a2[[i]], ttl=paste(a3[i],a4[i],sep = "_"), col_p=greenred(100),
              sep_vec =sep_vec)
  dev.off()
},a1=GT_GSVA_overlap_MEs,a2=ccp_cluster_assignments,a3="SKCM", a4=list(SKCM_clusters))

## association with AHR activation p.value=0.05 #####
TCGA_MEs_GSVA_GTs_pos_only <- map(TCGA_MEs_GSVA_GTs_not_null_pos_nodir, function(x){
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
## heatmaps
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
},a1=TCGA_modules_genes_pos_only_cpm_mats, a2=ccp_cluster_assignments,a3=tcga_names, a4=TCGA_clusters)


## Immune infiltration
cell_types_metagene_lists <- readRDS("/home/analyses/Projects_analyses/AHR_IL4I1/AHR/Publication_results/Resources/RDS/IPS_cell_types_gene_lists.rds")


## Generate the GSVA scores for the different immune populations
## estimate the GSVA scores for the different cell populations
GSVA_voom_Immune_cells <- map(list(v), function(a,b,c){
  glists <- map(b, function(x)x$metagene)
  gsva(expr = t(a$E),glists, parallel.sz=c,method="gsva")
},b=cell_types_metagene_lists,c=6)

saveRDS(GSVA_voom_Immune_cells,"../SKCM/SKCM_GSVA_voom_28_cell_types_IPS.rds")


## As per Charoentong et al 2017, after running a random forest classifier, the following cell types were
## identified to separate the patients with cytolytic activity from those with non-cytolytic activity.
## The active cells are active CD4+/CD8+ and Tem CD4/CD8, the supressors are MDSCs and T-regs.
## First plot heatmaps
## Load Cluster Assignments
## colors ####
cluster_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b",
                    "#e377c2","#7f7f7f","#bcbd22","#17becf")
heatmap_colors <- brewer.pal(5,"RdBu")[5:1]

####
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

GSVA_voom_Immune_cells_averaged_heatmaps <- map2(GSVA_voom_Immune_cells, list(clusters=ccp_cluster_assignments), function(a1,a2){
  df <- data.frame(cluster=a2, t(a1),  stringsAsFactors = F)
  colnames(df)[1] <- "clusters"
  split_levels <- levels(as.factor(df$cluster))
  split_dfs <- map(split_levels, function(b1){
    s_df <- df[df$cluster==b1,-1]
    colSums(s_df)
  })
  names(split_dfs) <- paste("G",split_levels, sep = "")
  res <- do.call(cbind,split_dfs)
})

heatmap_active_clusters <- c(1)

Averaged_GSVA_immune_heatmaps_all <- map(1, function(i,a1,a2,a3){
  annot_colors <- cluster_colors[1:ncol(a1[[i]])]
  names(annot_colors) <- colnames(a1[[i]])
  active_colors <- rep(cluster_colors[a3[i]], length(names(annot_colors)))
  names(active_colors) <- names(annot_colors)
  
  h1 <- HeatmapAnnotation(df = data.frame(cluster = names(annot_colors)),
                          col = list(cluster = annot_colors), show_legend = c(FALSE,FALSE))
  h2 <- HeatmapAnnotation(df = data.frame(Active = names(active_colors)),
                          col = list(Active = active_colors), show_legend = c(FALSE,FALSE))
  h3 <- Heatmap(a1[[i]], top_annotation = h1, bottom_annotation = h2, column_title = a2[i],
                  cluster_columns = FALSE, col = heatmap_colors, cluster_rows = FALSE,
                  row_names_side = "left",show_heatmap_legend = FALSE) 
  h3
},
a1=GSVA_voom_Immune_cells_averaged_heatmaps,
a2="SKCM",
a3=heatmap_active_clusters)

names(Averaged_GSVA_immune_heatmaps_all) <- "SKCM"

pdf("../SKCM/SKCM_immune_heatmap.pdf",width = 12, height = 8,pointsize = 18)
draw(Averaged_GSVA_immune_heatmaps_all$SKCM, gap=unit(6, "mm"), show_annotation_legend = TRUE)
dev.off()

###
IL4I1_voom <- map(list(v), function(x)x$E[,colnames(x$E)%in%c("IL4I1","IDO1","TDO2")])

IL4I1_GSVA_clusters_plots <- map(1, function(i,a1,a2,a3){
  df_test <- data.frame(gsva_sc=a1, clust=as.factor(a2[[i]]))
  comp_means_res <- compare_means(df_test, formula = gsva_sc~clust) %>%  as.data.frame()
  if(length(levels(df_test$clust))==2){
    comp_y_position <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
  } else if(length(levels(df_test$clust))==3){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
    comp_y_position <- seq(from = val_from,by = 0.2,length.out = 3)
  } else if(length(levels(df_test$clust))==4){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.3
    comp_y_position <- seq(from = val_from,by = 0.7,length.out = 6)
  }
  comp_means_res <-comp_means_res %>% mutate(.,y.position=comp_y_position) %>%
    mutate(sig=annt_fun(.$p.adj,length(.$p.adj)))
  p <- ggboxplot(df_test, x="clust", y="gsva_sc", color = "clust",
                 legend="right",palette = cluster_colors,ticks=FALSE,
                 scales="free_y",x.text.angle = 45,ylab = "log2 cpm",
                 title = a3[[i]])+
    stat_pvalue_manual(label = "p.adj",comp_means_res)+
    ggbeeswarm::geom_quasirandom(aes(color=clust))
  p
},
a1=IL4I1_voom[[1]][,1], a2=ccp_cluster_assignments, a3="IL4I1")

IDO1_GSVA_clusters_plots <- map(1, function(i,a1,a2,a3){
  df_test <- data.frame(gsva_sc=a1, clust=as.factor(a2[[i]]))
  comp_means_res <- compare_means(df_test, formula = gsva_sc~clust) %>%  as.data.frame()
  if(length(levels(df_test$clust))==2){
    comp_y_position <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
  } else if(length(levels(df_test$clust))==3){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
    comp_y_position <- seq(from = val_from,by = 0.2,length.out = 3)
  } else if(length(levels(df_test$clust))==4){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.3
    comp_y_position <- seq(from = val_from,by = 0.7,length.out = 6)
  }
  comp_means_res <-comp_means_res %>% mutate(.,y.position=comp_y_position) %>%
    mutate(sig=annt_fun(.$p.adj,length(.$p.adj)))
  p <- ggboxplot(df_test, x="clust", y="gsva_sc", color = "clust",
                 legend="right",palette = cluster_colors,ticks=FALSE,
                 scales="free_y",x.text.angle = 45,ylab = "log2 cpm",
                 title = a3[[i]])+
    stat_pvalue_manual(label = "p.adj",comp_means_res)+
    ggbeeswarm::geom_quasirandom(aes(color=clust))
  p
},
a1=IL4I1_voom[[1]][,2], a2=ccp_cluster_assignments, a3="IDO1")

TDO2_GSVA_clusters_plots <- map(1, function(i,a1,a2,a3){
  df_test <- data.frame(gsva_sc=a1, clust=as.factor(a2[[i]]))
  comp_means_res <- compare_means(df_test, formula = gsva_sc~clust) %>%  as.data.frame()
  if(length(levels(df_test$clust))==2){
    comp_y_position <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
  } else if(length(levels(df_test$clust))==3){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.1
    comp_y_position <- seq(from = val_from,by = 0.2,length.out = 3)
  } else if(length(levels(df_test$clust))==4){
    val_from <- max(df_test$gsva_sc,na.rm = TRUE)+0.3
    comp_y_position <- seq(from = val_from,by = 0.7,length.out = 6)
  }
  comp_means_res <-comp_means_res %>% mutate(.,y.position=comp_y_position) %>%
    mutate(sig=annt_fun(.$p.adj,length(.$p.adj)))
  p <- ggboxplot(df_test, x="clust", y="gsva_sc", color = "clust",
                 legend="right",palette = cluster_colors,ticks=FALSE,
                 scales="free_y",x.text.angle = 45,ylab = "log2 cpm",
                 title = a3[[i]])+
    stat_pvalue_manual(label = "p.adj",comp_means_res)+
    ggbeeswarm::geom_quasirandom(aes(color=clust))
  p
},
a1=IL4I1_voom[[1]][,3], a2=ccp_cluster_assignments, a3="TDO2")


IL4I1_voom_boxplots <- map2(IL4I1_voom, ccp_cluster_assignments, function(a1,a2){
  df <- data.frame(IL4I1=a1, cluster=paste("G",as.character(a2),sep = ""))
  })

df_melt <- IL4I1_voom_boxplots[[1]] %>% reshape::melt()

df_melt$variable <- gsub("IL4I1\\.","",df_melt$variable)



ggboxplot(df_melt, x="cluster", "value", facet.by = "variable", color = "cluster",
          palette = cluster_colors, scales="free")+ggbeeswarm::geom_quasirandom(aes(color=cluster))+stat_compare_means()

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
                 legend="right",palette = cluster_colors,ticks=FALSE,ylab = ttl,
                 scales="free_y",x.text.angle = 45
  )+stat_pvalue_manual(label = "p.adj", comp_means_res)+ggbeeswarm::geom_quasirandom(aes(color=cluster))
  p
}

IL4I1_boxplots_list <- map2(IL4I1_voom_boxplots,names(TCGA_voom),dfs_to_plot_FUN)
