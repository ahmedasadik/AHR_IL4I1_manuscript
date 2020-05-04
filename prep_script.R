#!/usr/bin/env Rscript

## installation of requirements
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
}

ipak_bc <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    BiocManager::install(new.pkg)
}

dep_packs <- c("BiocManager","limma", "purrr", "oligo", "pd.hugene.2.0.st", "RColorBrewer",
               "affycoretools", "biomaRt", "GSEABase", "clusterProfiler", "ggplot2", "gplots",
               "org.Hs.eg.db", "sva", "parallel", "ggpubr", "reshape", "topGO", "remotes",
               "WGCNA", "GEOquery", "edgeR", "rJava", "circlize", "ComplexHeatmap", "parallel",
               "dplyr", "xlsx", "corrplot", "psych"
               )

ipak(dep_packs)
ipak_bc(dep_packs)

remotes::install_github("ahmedasadik/AffyGEx", auth_token="5d53366564bec63272a087be6228eb7286085d8a")
