#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript/AHR_signature/
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This files describes how the AHR signature was generated
################################################

## Load libraries
library(purrr)
library(biomaRt)
library(parallel)

## Source the functions and parameters files
source("../functions_and_parameters.R")

#####################
## NLP AHR targets ##
#####################
## These are the targets that were generated as a result of the text mining approach
nlp_AHR <- read.delim("../Resources/ahr_targets_filtered_abundances.csv")
colnames(nlp_AHR) <- c("gene", "freq")

## Extracting unique AHR targets and converting them to upper case
nlp_AHR_genes <- nlp_AHR$gene %>% as.character() %>% toupper() %>% unique()

## Renaming the genes that require their name changed
nlp_AHR_genes[c(51,80,122,207,423,447)] <- c("AR", "TH","F3", "HP","INFA","JUN")

#######################
## Ahmed AHR targets ##
#######################
## These are the lists of differentially regulated genes acquired from curated datasets
# GEO datasets
sheet_path <- "../Resources/selected_AHR_studies.xlsx"
datasheets <- xlsx::loadWorkbook(sheet_path)
sheets <- xlsx::getSheets(datasheets)

ahr_signalsheets <- map2(sheet_path, names(sheets), function(a,b){read.xlsx(file = a, sheetName = b, stringsAsFactors=F)})
names(ahr_signalsheets) <- names(sheets)

genes <- map(ahr_signalsheets, function(x){x[,1]}) %>% unlist() %>% unique()

all_signals <- map(genes, function(gene_sym){
  idcs <- map(ahr_signalsheets, function(x){
    idx <-which(gene_sym == x$Gene)
    if(length(idx)!=0){
      return(idx)
    }else{
      return(NA)
    }
  }) %>% unlist()
  to_df <- as.data.frame(matrix(nrow = 1, ncol = 4))
  for(i in seq_along(idcs)){
    if(!is.na(idcs[i])){
      to_df[i,] <- ahr_signalsheets[[i]][idcs[i],]
    } else {
      to_df[i,] <- NA
    }
  }
  to_df <- to_df[complete.cases(to_df)==TRUE,]
  gen_state <-c(state=paste(to_df[,2], collapse = ";"))
  gen_conditions <- c(condition=paste(to_df[,3], collapse = ";"))
  gen_pmid <- c(pmid=paste(to_df[,4], collapse = ";"))
  final_df <- data.frame(gene_sym, gen_state, gen_conditions,gen_pmid, stringsAsFactors = F)
})

signal_df <- do.call("rbind", all_signals)

write.table(signal_df, "../Results/AHR_signature/signal_df.txt", sep = "\t")

## Get the gene symbols of the features annotated with ensemble transcript IDs
signal_df_genes_enst <- signal_df$gene_sym[grep("^ENST", signal_df$gene_sym)]

listMarts()
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters_ensembl <- listFilters(ensembl)
attributes_ensembl <- listAttributes(ensembl)
signal_df_genes_enst_genes <- getBM(attributes=c('ensembl_transcript_id', 'entrezgene', 'hgnc_symbol'),
                                    filter='ensembl_transcript_id', values = signal_df_genes_enst,
                                    mart = ensembl)
signal_df_genes_all <- signal_df$gene_sym %>% .[-grep("^ENST", .)] %>% c(.,signal_df_genes_enst_genes$hgnc_symbol[!is.na(signal_df_genes_enst_genes$entrezgene)])

############
## Sygnal ##
############
## These are the targets that have AHR binding sites as per the Sygnal database (tfbsdb.systemsbiology.net)
# retrieve the AHR motif targets
Sygnal_AHR <- list.files("../Resources/Sygnal/", full.names = T)
Sygnal_AHR_files <- map(Sygnal_AHR, read.delim, sep=",", stringsAsFactors=F)
Sygnal_AHR_eIDs <- Sygnal_AHR_files %>% map(.,function(x){x[,1]}) %>% unlist() %>% as.character()

# get the gene symbols of the AHR target genes
Sygnal_AHR_genes <- unlist(mget(unique(Sygnal_AHR_eIDs), org.Hs.egSYMBOL,ifnotfound=NA))

###################################
## Annotation update as per HGNC ##
###################################
# This was downloaded directly from HGNC on 04-12-2017
hgnc <- read.delim("../Resources/human_hgnc_annotation_file.txt")

## This is performed for all three target types
cl <- makeCluster(no_cores)
clusterExport(cl,varlist = c("hgnc"))

## Signal genes
signal_df_genes_all_hgnc_annot <- t(parSapply(cl, signal_df_genes_all, FUN = annot_query_FUN,
                                              sym = hgnc$symbol, eid = hgnc$entrez_id, syn = hgnc$alias_symbol,
                                              psym = hgnc$prev_symbol))
signal_df_genes_all_hgnc_annot <- signal_df_genes_all_hgnc_annot[-c(which(duplicated(rownames(signal_df_genes_all_hgnc_annot))==TRUE)),]

## Sygnal genes
Sygnal_AHR_genes_hgnc_annot <- t(parSapply(cl, Sygnal_AHR_genes, FUN = annot_query_FUN,
                                              sym = hgnc$symbol, eid = hgnc$entrez_id, syn = hgnc$alias_symbol,
                                              psym = hgnc$prev_symbol))
## nlp genes
nlp_AHR_genes_hgnc_annot <- t(parSapply(cl, nlp_AHR_genes, FUN = annot_query_FUN,
                                           sym = hgnc$symbol, eid = hgnc$entrez_id, syn = hgnc$alias_symbol,
                                           psym = hgnc$prev_symbol))

stopCluster(cl)

### Signal genes proofing after HGNC 
## signal_df_genes_all_hgnc_annot
signal_df_genes_all_hgnc_annot_na <- as.data.frame(signal_df_genes_all_hgnc_annot) %>% mutate_all(funs(empty_as_na))
rownames(signal_df_genes_all_hgnc_annot_na) <- rownames(signal_df_genes_all_hgnc_annot)
problem_rows_signal <- rownames(signal_df_genes_all_hgnc_annot_na)[which(is.na(signal_df_genes_all_hgnc_annot_na[,1]))]

## retrieve the original entries from the signal_df dataframe and check the PMIDs
problem_df_signal <- signal_df[signal_df$gene_sym %in% problem_rows_signal,]
signal_to_correct <- c("FBX027", "DDlT4", "SLC3SG5", "C60rf132","TNAP", "TMDCII", "SPANXB2", "PLXNA4A", "PLXN3",
                       "NALP1", "C9orf79", "C9orf14", "C21orf42", "LOC285593", "C14orf166B", "SFRS2IP",
                       "HDGFRP3", "Mar-10", "Sep-10", "HAS2AS", "C21orf130", "C21orf15", "C4orf12",
                       "C21orf29", "KLRA1")

## correct the entries
signal_correct_gs <- c("FBXO27", "DDIT4", "SLC35G5", "C6orf132", "ALPL","ADAM5", "SPANXB1", "PLXNA4", "PLXNA3",
                       "NLRP1", "SPATA31E1", "LINC00032", "LINC00158", "LOC285593", "LRRC74A", "SCAF11",
                       "HDGFL3", "MARCH10", "SEPT10", "HAS2-AS1", "LINC00323", "CYP4F29P", "WDFY3-AS2",
                       "TSPEAR", "KLRA1P")

signal_correct_EIDS <- c("126433", "54541", "83650", "647024", "249", "255926", "728695", "91584", "55558",
                         "22861", "286234", "158635", "54072", "285593", "145497", "9169", "50810", "162333",
                         "151011", "594842", "284835", "54055", "404201", "54084", "10748")

signal_df_genes_all_hgnc_annot_na$V1[(rownames(signal_df_genes_all_hgnc_annot_na)%in%signal_to_correct)] <- signal_correct_gs
signal_df_genes_all_hgnc_annot_na$V2[(rownames(signal_df_genes_all_hgnc_annot_na)%in%signal_to_correct)] <- signal_correct_EIDS

MLL4_idx <- which(rownames(signal_df_genes_all_hgnc_annot_na)=="MLL4")
MLL4_gs <- strsplit(signal_df_genes_all_hgnc_annot_na[MLL4_idx,1],";")%>% unlist()
MLL4_eids <- strsplit(signal_df_genes_all_hgnc_annot_na[MLL4_idx,2],";")%>% unlist()

signal_df_genes_all_hgnc_annot_na <- signal_df_genes_all_hgnc_annot_na[-MLL4_idx,]
signal_df_genes_all_hgnc_annot_na <- rbind(signal_df_genes_all_hgnc_annot_na, data.frame(V1=MLL4_gs, V2=MLL4_eids))
signal_df_genes_all_hgnc_annot_na <- signal_df_genes_all_hgnc_annot_na[!is.na(signal_df_genes_all_hgnc_annot_na$V1),]
signal_df_genes_all_hgnc_annot_na <- signal_df_genes_all_hgnc_annot_na[!is.na(signal_df_genes_all_hgnc_annot_na$V2),]
signal_df_genes_all_hgnc_annot_na <- signal_df_genes_all_hgnc_annot_na[which(duplicated(signal_df_genes_all_hgnc_annot_na$V1)==FALSE),]
write.table(signal_df_genes_all_hgnc_annot_na,"../Results/AHR_signature/AHR_genes_datasets.txt", sep = "\t")

### Sygnal genes proofing after HGNC
## Sygnal_AHR_genes_hgnc_annot
Sygnal_AHR_genes_hgnc_annot_na <- as.data.frame(Sygnal_AHR_genes_hgnc_annot) %>% mutate_all(funs(empty_as_na))
rownames(Sygnal_AHR_genes_hgnc_annot_na) <- rownames(Sygnal_AHR_genes_hgnc_annot)
# find rows with no gene names
problem_rows_Sygnal <- rownames(Sygnal_AHR_genes_hgnc_annot_na)[which(is.na(Sygnal_AHR_genes_hgnc_annot_na[,1]))]
# correct the values of these values
Sygnal_AHR_genes_hgnc_annot_na[problem_rows_Sygnal,2] <- problem_rows_Sygnal
Sygnal_AHR_genes_hgnc_annot_na[problem_rows_Sygnal,1] <- c("LOC100287036")
write.table(Sygnal_AHR_genes_hgnc_annot_na, "../Results/AHR_signature/AHR_genes_Sygnal.txt", sep = "\t")

### nlp genes proofing after HGNC
## nlp_AHR_genes_hgnc_annot
nlp_AHR_genes_hgnc_annot_na <- as.data.frame(nlp_AHR_genes_hgnc_annot) %>% mutate_all(funs(empty_as_na))
# find rows with no gene names
problem_rows_nlp <- rownames(nlp_AHR_genes_hgnc_annot)[which(is.na(nlp_AHR_genes_hgnc_annot_na[,1]))]

# proofing gene entries and removing those with no human orthologous
nlp_ok_idx <- which(rownames(nlp_AHR_genes_hgnc_annot)=="C21ORF33" |
                         rownames(nlp_AHR_genes_hgnc_annot)=="LOC694944" |
                         rownames(nlp_AHR_genes_hgnc_annot)=="MAMU-DMA" |
                         rownames(nlp_AHR_genes_hgnc_annot)=="UGT1A2" |
                         rownames(nlp_AHR_genes_hgnc_annot)=="ANC" |
                         rownames(nlp_AHR_genes_hgnc_annot)=="CPX" |
                         rownames(nlp_AHR_genes_hgnc_annot)=="CYP1D1" |
                         rownames(nlp_AHR_genes_hgnc_annot)=="DBA2")

nlp_AHR_genes_hgnc_annot_na[nlp_ok_idx,1] <- c("C21orf33","SAA1","HLA-DMA","UGT1A3","ANC", "TBX22","CYP1D1P", "DBA2")
nlp_AHR_genes_hgnc_annot_na[nlp_ok_idx,2] <- c("8209", "6288", "3108", "54659", "8066", "50945", "100133307", "114086")
nlp_AHR_genes_hgnc_annot_na2 <- nlp_AHR_genes_hgnc_annot_na[!is.na(nlp_AHR_genes_hgnc_annot_na$V1),]

DEPP_idx <- grep("DEPP", nlp_AHR_genes_hgnc_annot_na2$V1)
MT1E_idx <- grep("MT1E", nlp_AHR_genes_hgnc_annot_na2$V1)

DEPP_gs <- strsplit(nlp_AHR_genes_hgnc_annot_na2$V1[DEPP_idx],";")%>% unlist()
DEPP_eids <- strsplit(nlp_AHR_genes_hgnc_annot_na2$V2[DEPP_idx],";")%>% unlist()

MT1E_gs <- strsplit(nlp_AHR_genes_hgnc_annot_na2$V1[MT1E_idx],";")%>% unlist()
MT1E_eids <- strsplit(nlp_AHR_genes_hgnc_annot_na2$V2[MT1E_idx],";")%>% unlist()

nlp_AHR_genes_hgnc_annot_na2 <- nlp_AHR_genes_hgnc_annot_na2[-c(DEPP_idx, MT1E_idx),]
nlp_AHR_genes_hgnc_annot_na2 <- rbind(nlp_AHR_genes_hgnc_annot_na2, data.frame(V1=DEPP_gs,V2=DEPP_eids),
                                         data.frame(V1=MT1E_gs, V2=MT1E_eids), append=T)

write.table(nlp_AHR_genes_hgnc_annot_na2, "../Results/AHR_signature/AHR_genes_nlp.txt", sep = "\t")

## Combining signal and Sygnal dataframes
AHR_genes_datasets <- c(Sygnal_AHR_genes_hgnc_annot_na$V1, signal_df_genes_all_hgnc_annot_na$V1)%>% .[which(duplicated(.)==FALSE)]

AHR_genes_datasets <- dplyr::full_join(Sygnal_AHR_genes_hgnc_annot_na,signal_df_genes_all_hgnc_annot_na)

write.table(as.data.frame(AHR_genes_datasets), "../Results/AHR_signature/AHR_genes_ds_S.txt", sep = "\t", col.names = F)

## Overlapping genes to be used as a signature for AHR
overlapping_genes <- AHR_genes_datasets[AHR_genes_datasets$V1 %in% nlp_AHR_genes_hgnc_annot_na2$V1,]
colnames(overlapping_genes) <- c("Gene", "EID")
write.table(overlapping_genes, "../Results/AHR_signature/overlapping_AHR_signature_genes.txt", sep = "\t")
