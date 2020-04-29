# This script describes updating the annotations of the metagenes used for identifying infiltrtaing cells
# as per the Charoentong et al 2017 Cell Reports paper.

## Load libraries
library(purrr)
library(parallel)
library(dplyr)

## Functions
# Change blanks to NA
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

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

## HGNC
hgnc <- read.delim("./Resources/HGNC/hgnc_complete_set_04122017.txt", stringsAsFactors = F)

## Read the immune metagenes
metagenes <- read.csv("./Tables/IPS_cell_metagenes.csv", stringsAsFactors = F)
colnames(metagenes) <- metagenes[2,]
metagenes <- metagenes[-c(1:2),]

## Update the annotations
no_cores <- detectCores() - 2
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

saveRDS(cell_types_metagene_lists, "./RDS/IPS_cell_types_gene_lists.rds")
