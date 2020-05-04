#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript/AHR_signature/
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes updating the annotations of the metagenes used for identifying infiltrtaing cells
## as per the Charoentong et al 2017 Cell Reports paper (the metagenes were downloaded from the publication supplements)
################################################

## Source the functions and parameters files
source("./functions_and_parameters.R")

## Load libraries
library(purrr)
library(parallel)
library(dplyr)

## HGNC
hgnc <- read.delim("./Resources/human_hgnc_annotation_file.txt", stringsAsFactors = F)

## Read the immune metagenes
metagenes <- read.csv("./Resources/IPS_cell_metagenes.csv", stringsAsFactors = F)
colnames(metagenes) <- metagenes[2,]
metagenes <- metagenes[-c(1:2),]

## Update the annotations
cl <- makeCluster(no_cores)
clusterExport(cl,varlist = c("hgnc"))
metagene_hgnc_annot <- t(parSapply(cl, metagenes$Metagene, FUN = annot_query_FUN, 
                                   sym = hgnc$symbol, eid = hgnc$entrez_id, syn = hgnc$alias_symbol, 
                                   psym = hgnc$prev_symbol))
stopCluster(cl)

## change blanks into NA
metagene_hgnc_annot_na <- as.data.frame(metagene_hgnc_annot) %>% mutate_all(funs(empty_as_na))
rownames(metagene_hgnc_annot_na) <- rownames(metagene_hgnc_annot)

## Update gene symbols to the recent release (these genes were not found in the HGNC release used)
metagene_hgnc_annot_na[which(rownames(metagene_hgnc_annot_na)=="HDGFRP2"),1] <- "HDGFL2"
metagene_hgnc_annot_na[which(rownames(metagene_hgnc_annot_na)=="HDGFRP2"),2] <- "84717"
metagene_hgnc_annot_na[which(rownames(metagene_hgnc_annot_na)=="ClQA"),1] <- "C1QA"
metagene_hgnc_annot_na[which(rownames(metagene_hgnc_annot_na)=="ClQA"),2] <- "712"
metagene_hgnc_annot_na[which(rownames(metagene_hgnc_annot_na)=="ClQB"),1] <- "C1QB"
metagene_hgnc_annot_na[which(rownames(metagene_hgnc_annot_na)=="ClQB"),2] <- "713"

metagene_hgnc_annot_na <- data.frame(metagene=metagene_hgnc_annot_na$V1,EID=metagene_hgnc_annot_na$V2,
                                     cell_type=metagenes$`Cell type`,Immunity_type=metagenes$Immunity)

metagene_hgnc_annot_na <- metagene_hgnc_annot_na[which(duplicated(metagene_hgnc_annot_na$metagene)==FALSE),]
metagene_hgnc_annot_na$cell_type <- gsub(" ","_",metagene_hgnc_annot_na$cell_type)
cell_types <- levels(as.factor(metagene_hgnc_annot_na$cell_type))

## Generate a list of the updated anotated gene signatures for the each cell type
cell_types_metagene_lists <- map(cell_types, function(a){metagene_hgnc_annot_na[metagene_hgnc_annot_na$cell_type==a,]})
names(cell_types_metagene_lists) <- cell_types

saveRDS(cell_types_metagene_lists, "./Results/RDS/IPS_cell_types_gene_lists.rds")
