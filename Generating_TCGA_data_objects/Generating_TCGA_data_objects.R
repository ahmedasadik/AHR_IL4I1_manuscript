#!/usr/bin/env Rscript

#################################################
## Project: AHR LAAO
## Origin: https://github.com/ahmedasadik/Project_AHR_LAAO/tree/master/AHR_scripts
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This file describe the TCGA data objects that were used throughout the project.
## This includes downloading and processing the data
## Only harmonized TCGA RNA-seq data was used mapped to hg38
## Both HTseq_counts and FPKMs were downloaded
## Only primary tumors were retained and all others (including metastasis) were filtered out
## FPKMs were further converted to TPMs
#################################################

## Load libraries
library(parallel)
library(TCGAbiolinks)
library(purrr)
library(edgeR)
library(miceadds)

# Functions imported from functions_and_parameters.R
source("../functions_and_parameters.R")

###############################
## Downloading the TCGA data ##
###############################
#######################
## HTSeq-Counts data ##
#######################
## This downloads the counts RData objects into the TCGA_counts/RData folder
setwd("../TCGAbiolinks_downloads/TCGA_counts")
cl <- makeCluster(no_cores)
parLapply(cl=cl, project_ids, counts_download_FUN)
stopCluster(cl)
setwd("../")

## Preparing the data into DGElists 
# Extracting the assay data and creating DGElists
# The counts data were filtered for only primary tumors

tcga_datasests_paths <- list.files("../TCGAbiolinks_downloads/TCGA_counts/RDats/", full.names = TRUE)
tcga_datasests_names <- list.files("../TCGAbiolinks_downloads/TCGA_counts/RDats/") %>% gsub("Tr.*","",.) %>% gsub("\\-", "_",.)

# read all tcga summarized experiment data
TCGA_datasest <- map(tcga_datasests_paths, load.Rdata)
names(TCGA_datasest) <- tcga_datasests_names

# generate the count data
TCGA_cntdata <- map(TCGA_datasest, SummarizedExperiment::assay)

# generate the coldata
TCGA_coldata <- map(TCGA_datasest, SummarizedExperiment::colData)
names(TCGA_coldata) <- names(TCGA_cntdata)

# generate the rowdata (this is identical for all samples across all tumors)
TCGA_rowdata <- map(TCGA_datasest, SummarizedExperiment::rowData)
names(TCGA_rowdata) <- names(TCGA_cntdata)

# generate the TP sample indeces
samplesTP_idcs <- map(TCGA_coldata, function(cnts){TCGAbiolinks::TCGAquery_SampleTypes(barcode = cnts$barcode,typesample = c("TP"))}) 
names(samplesTP_idcs) <- names(TCGA_coldata)

# filter TCGA counts and coldata for TP samples only
tcga_cnts_TP <- map2(TCGA_cntdata,samplesTP_idcs,function(x,y){x[,match(y, colnames(x))]})
tcga_coldata_TP <- map2(TCGA_coldata,samplesTP_idcs,function(x,y){x[match(y, x$barcode),]})

# Creating DGE objects for all tumors
TCGA_DGE <- map2(tcga_cnts_TP, tcga_coldata_TP, function(x,y,z){DGEList(counts = x, samples = y, genes = z)}, z=TCGA_rowdata[[1]])
saveRDS(TCGA_DGE, "../Results/RDS/TCGA_DGE.rds")

## Generating an object only for TCGA counts data that were converted to log2(cpm+1)
TCGA_counts <- map(TCGA_datasest, function(x){
  cnt_data <-  SummarizedExperiment::assay(x)
  samplesTP <- TCGAbiolinks::TCGAquery_SampleTypes(barcode = colnames(cnt_data),typesample = c("TP")) 
  counts_TP <- cnt_data[,match(samplesTP, colnames(cnt_data))]
  row_data <- SummarizedExperiment::rowData(x)
  rownames(counts_TP) <- row_data$external_gene_name[match(rownames(counts_TP), row_data$ensembl_gene_id)]
  counts_cpm <- edgeR::cpm(counts_TP, prior.count=1, log=T)
})
saveRDS(TCGA_counts, "../Results/RDS/TCGA_counts.rds")

## Generating an object only for TCGA coldata data that were converted to log2(cpm+1)
TCGA_colData <- map2(TCGA_counts, TCGA_datasest, function(cnts_mat, col_df){
  col_data <-  SummarizedExperiment::colData(col_df)
  new_ColData <- col_data[match(colnames(cnts_mat), col_data$barcode),]
})
saveRDS(TCGA_colData, "../Results/RDS/Results/TCGA_coldata.rds")

## Generating a normalized counts DGElist using TMM normalization and voom. The data was filtered first.
# Voom to prepare the counts for WGCNA
TCGA_DGE_voom <- map(TCGA_DGE, function(dge){
  keep <- filterByExpr(dge)
  dge <- dge[keep,,keep.lib.size=FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, plot=FALSE)
  v
})

# Saving the TCGA voomed DGE to disk
saveRDS(TCGA_DGE_voom, "../Results/RDS/TCGA_DGE_voom.rds")

## Generate the voomed TCGA dataset that will be used for WGCNA
TCGA_DGE_voom_annots <- map(TCGA_DGE_voom, function(dge){
  rownames(dge$E) <- dge$genes$external_gene_name [match(rownames(dge$E), dge$genes$ensembl_gene_id)]
  dge$E <- t(dge$E)
  dge
})
saveRDS(TCGA_DGE_voom_annots, "../Results/RDS/TCGA_DGE_voom_annot.rds")

#####################
## HTSeq-FPKM data ##
#####################
## This downloads the FPKM files in the TCGA_FPKM/RData that will be read later to generate the TPM Object
setwd("../TCGAbiolinks_downloads/TCGA_FPKM")
cl <- makeCluster(no_cores)
parLapply(cl=cl, project_ids, counts_download_FUN)
stopCluster(cl)
setwd("../")
# List files of TCGA datasets
FPKM_tcga_datasets_paths <- list.files("../TCGAbiolinks_downloads/TCGA_FPKM/RDats", full.names = T)
FPKM_tcga_datasets_names <- list.files("../TCGAbiolinks_downloads/TCGA_FPKM/RDats") %>% gsub("Tr.*","",.) %>% gsub("\\-", "_",.)

# Loading all HT-Seq FPKMs datasets
FPKM_TCGA_datasets <- map2(FPKM_tcga_datasets_paths, FPKM_tcga_datasets_names, load.Rdata)
names(FPKM_TCGA_datasets) <- FPKM_tcga_datasests_names

# generate the FPKM data
TCGA_fpkmdata <- map(FPKM_TCGA_datasest, SummarizedExperiment::assay)

# generate the coldata
FPKM_TCGA_coldata <- map(FPKM_TCGA_datasest, SummarizedExperiment::colData)
names(FPKM_TCGA_coldata) <- names(TCGA_fpkmdata)

# generate the rowdata (this is identical for all samples across all tumors)
FPKM_TCGA_rowdata <- map(FPKM_TCGA_datasest, SummarizedExperiment::rowData)
names(FPKM_TCGA_rowdata) <- names(TCGA_fpkmdata)

# generate the TP sample indeces
FPKM_samplesTP_idcs <- map(FPKM_TCGA_coldata, function(cnts){TCGAbiolinks::TCGAquery_SampleTypes(barcode = cnts$barcode, typesample = c("TP"))}) 
names(FPKM_samplesTP_idcs) <- names(FPKM_TCGA_coldata)

# filter TCGA FPKMs and coldata for TP samples only
tcga_fpkm_TP <- map2(TCGA_fpkmdata, FPKM_samplesTP_idcs, function(x,y){x[,match(y, colnames(x))]})
FPKM_tcga_coldata_TP <- map2(FPKM_TCGA_coldata, FPKM_samplesTP_idcs,function(x,y){x[match(y, x$barcode),]})

# Creating DGE objects for all tumors
TCGA_fpkm_DGE <- map2(tcga_fpkm_TP, FPKM_tcga_coldata_TP, function(x,y,z){DGEList(counts = x, samples = y, genes = z)}, z=FPKM_TCGA_rowdata[[1]])
saveRDS(TCGA_fpkm_DGE, "../Results/RDS/TCGA_fpkms_all.rds")

##############
## TPM data ##
##############
## Converting FPKMs to TPMs, which can be used later for Cibersort or other plots if need be.
TCGA_TPMs <- map(TCGA_fpkms_DGE, safely(function(dge){
  dge$counts <- fpkmToTpm(dge$counts)
  colnames(dge$counts) <- paste("Pt", 1:length(colnames(dge$counts)))
  genes <- dge$genes$external_gene_name[match(rownames(dge$counts), dge$genes$ensembl_gene_id)]
  tpm_df <- data.frame(GeneSymbol=genes, dge$counts, stringsAsFactors = F)
}))
saveRDS(TCGA_TPMs, "../Results/RDS/TCGA_TPMs.rds")
