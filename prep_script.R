## installation of requirements

dir.create("./Results")
dir.create("./Results/AHR_signature")
dir.create("./Results/AHR_signature_validations")
dir.create("./Results/AHR_Trp_metabolism_associations")
dir.create("./Results/")

library(purrr)
library(biomaRt)
library(parallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
library(dplyr)