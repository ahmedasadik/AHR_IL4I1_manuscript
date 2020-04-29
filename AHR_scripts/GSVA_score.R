# This file explains the generation of a GSVA score for the AHR signature

## Load libraries
library(edgeR)
library(GSVA)

# Loading TCGA DGE object raw counts before normalization as DGElists
TCGA_DGE <- readRDS("./RDS/TCGA/TCGA_DGE.rds")

# Read the AHR signature file
overlapping_genes <- read.delim("./Signature/overlapping_AHR_signature.txt", sep = "\t", stringsAsFactors = F)

# Extracting the AHR ENSG_IDs
AHR_ensg <- TCGA_DGE$TCGA_ACC$genes[TCGA_DGE$TCGA_ACC$genes$external_gene_name %in% overlapping_genes$Gene,]

# generate the GSVA enrichment scores of all TCGA tumors
TCGA_gsva <- map(TCGA_DGE, safely(function(dge){(gsva(cpm(dge$counts, prior.count = 1, log = T),list(AHR_ensg$ensembl_gene_id), parallel.sz=6))}))
saveRDS(TCGA_gsva, "./RDS/TCGA_GSVA_scores_safely.rds")
