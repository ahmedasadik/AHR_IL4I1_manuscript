#!/usr/bin/env Rscript

## script to automate download of required libraries/packages/dependencies
## the CRAN mirror is set to http://cran.us.r-project.org, which can be changed
## add ", INSTALL_opts = '--no-lock'" if there are problems installing libraries

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.us.r-project.org')
}

ipak_bc <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg)
}

dep_packs <- c("BiocManager","limma", "purrr", "oligo", "pd.hugene.2.0.st", "RColorBrewer", "stringr","illuminaHumanv4.db","illuminaHumanv3.db",
               "affycoretools", "biomaRt", "GSEABase", "clusterProfiler", "ggplot2", "gplots", "hgu133plus2.db", "gcrma",
               "org.Hs.eg.db", "sva", "parallel", "ggpubr", "reshape", "topGO", "remotes", "lumi",
               "WGCNA", "GEOquery", "edgeR", "rJava", "circlize", "ComplexHeatmap", "parallel",
               "dplyr", "xlsx", "corrplot", "psych", "Hmisc", "factoextra", "FactoMineR",
               "huex10stprobeset.db", "R.utils", "extrafonts", "ggbeeswarm", "survival", "survminer",
               "GSVA", "globaltest", "TCGAbiolinks", "miceadds", "fmsb", "SummarizedExperiment",
               "pheatmap", "MultiAssayExperiment", "MOFA", "gridExtra", "grid"
               )

ipak(dep_packs)
ipak_bc(dep_packs)

remotes::install_github("ahmedasadik/AffyGEx", auth_token="5d53366564bec63272a087be6228eb7286085d8a", ref="v1.0.0")
