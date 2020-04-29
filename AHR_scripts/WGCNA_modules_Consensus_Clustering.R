# This script describes performing Consensus Clustering of the modules associating with AHR activation using Kmeans

## We need to cluster the patients into groups based on their levels of AHR activation. Deciding on the
## right number of clusters is important as well as the selection of the clustering algorithm. We will use
## as features the MEs associating with AHR activation (positive and negative associations). Defining the
## right number of clusters is performed by assessing the CDF plots, track plots, consensus heatmaps,
## cluster identity and the number of clusters using the elbow method.

## Load libraries
library(purrr)
library(ConsensusClusterPlus)
library(parallel)

# Module Eigen genes
GT_GSVA_overlap_MEs <- readRDS("./RDS/GT_GSVA_overlap_MEs.rds")

## CCP_fun 1000 subsamplings kmeans and distance 1-pearson
ccp_fun_km_dist <- function(mods, dirs){
  d_m <- 1-cor(as.matrix(mods)) %>% as.dist()
  ccp_run <- ConsensusClusterPlus(d_m, maxK = 10, reps = 1000,
                                  clusterAlg = "kmdist", writeTable = FALSE, seed = 0860,
                                  title = paste(dirs,"_ccp", sep = ""), plot="pdf")
  calcIC_run <- calcICL(ccp_run,plot = "pdf", title = paste(dirs,"_calcIC", sep = ""))
  res <- list(ccp=ccp_run, calcIC=calcIC_run)
}

# CCP run
no_cores <- detectCores() -1
cl <- makeCluster(no_cores)
clusterEvalQ(cl, {library(ConsensusClusterPlus);library(purrr)})
clusterExport(cl,varlist = c("ccp_fun_km_dist"))
idcs_list <- lapply(1:32, function(x)x)
GT_GSVA_overlap_MEs_t <- map(GT_GSVA_overlap_MEs, t)
ccp_pos_dist <- parLapply(cl, idcs_list, function(x, a1, a2){
  ccp_fun_km_dist(a1[[x]],a2[[x]])
}, a1=GT_GSVA_overlap_MEs_t, a2=paste("GTs", names(GT_GSVA_overlap_MEs_t),"dist", sep = "_"))
gc()
stopCluster(cl)

# Save
saveRDS(ccp_pos_dist, "./RDS/ccp_pos_dist.rds")
