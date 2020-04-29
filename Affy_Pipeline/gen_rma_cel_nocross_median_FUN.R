## This function generates a normalized and filtered rma object
gen_rma_cel_nocross_median_FUN <- function (raw_set, sd_val = 0, snames, ct_off){
  rma_cel <- rma(raw_set, target = "probeset")
  cat("Adding feature data at probeset level", fill = T)
  featureData(rma_cel) <- getNetAffx(rma_cel, "probeset")
  split_gene_assignment_FUN <- function(gene_assignment) {
    genes_split <- strsplit(gene_assignment, "\\//")
    genes_split[[1]][2] <- gsub(" ", "", genes_split[[1]][2])
    genes_split[[1]][2]
  }
  split_gene_assignment_ENST_FUN <- function(mrna_assignment) {
    str_idx <- grep("ENST", strsplit(mrna_assignment, "\\//")[[1]])
    str_split <- strsplit(mrna_assignment, "\\//")
    str_split[[1]][str_idx]
    sub_enst <- gsub("^\\/ ", "", str_split[[1]][str_idx])
    paste(gsub(" $", "", sub_enst), collapse = ";")
  }
  cat("Creating negative control index (used for all filtering calculations)", fill = T)
  neg_idx <- which(rma_cel@featureData@data$probesettype == "normgene->intron")
  filt_rule_sd_FUN <- function(x, n, l = 0) {
    if (l == 0) {
      neg_mean <- mean(exprs(rma_cel[neg_idx, n]))
      fil_gene_names <- names(x[x >= neg_mean])
    } else if (l > 0) {
      neg_mean <- mean(exprs(rma_cel[neg_idx, n]))
      neg_sd <- sd(exprs(rma_cel[neg_idx, n]))
      neg_filt <- neg_mean + l * neg_sd
      fil_gene_names <- names(x[x >= neg_filt])
    }
  }
  cat(paste0("Retaining probesets that are equal to or greater than ", 
             sd_val, " sd + the mean of negative controls"), fill = T)
  cat("Creating expression set object", fill = T)
  rma_exprs <- exprs(rma_cel)
  n_sam <- length(rma_cel@phenoData@data$sample_id)
  n <- list()
  for (i in 1:n_sam) {
    n[[i]] <- i
  }
  args_list <- list()
  for (i in 1:n_sam) {
    args_list[[i]] <- list(x = rma_exprs[, i], n = i, l = sd_val)
  }
  cat("Matching unique probesets to the different chips", fill = T)
  filt_prbsets_2sd <- invoke_map(list(filt_rule_sd_FUN), args_list)
  prbsets_2sd <- unlist(filt_prbsets_2sd)
  unique_idx_2sd <- unique(prbsets_2sd)
  unique_prbsets_2sd <- prbsets_2sd[match(unique_idx_2sd, prbsets_2sd)]
  cat("Matching filtered probesets", fill = T)
  filt_prbset_match_2sd <- sapply(filt_prbsets_2sd, function(x) {match(unique_prbsets_2sd, x)})
  cat(paste0("Filtering out probesets if present in less than ", 
             c(100 - ct_off * 100), "% of chips"), fill = T)
  filt_df_comp_2sd <- as.data.frame(matrix(nrow = length(unique_prbsets_2sd), ncol = 1 + n_sam))
  filt_df_comp_2sd[, 1] <- unique_prbsets_2sd
  filt_df_comp_2sd[, 2:c(1 + n_sam)] <- filt_prbset_match_2sd
  colnames(filt_df_comp_2sd) <- c("probesetid", snames)
  filt_df_comp_2sd$noNA <- unlist(apply(filt_df_comp_2sd, 1, function(x) {
    length(which(is.na(x) == TRUE))}))
  to_use_prbsetid_2sd <- filt_df_comp_2sd$probesetid[which(filt_df_comp_2sd$noNA <= round(ct_off * c(length(snames))))]
  idx_2sd <- match(to_use_prbsetid_2sd, rma_cel@featureData@data$probesetid)
  rma_cel_2sd <- rma_cel[idx_2sd]
  cat("Removing all remaining control or unannotated probesets", fill = T)
  idx_main_2sd <- which(rma_cel_2sd@featureData@data$probesettype == "main")
  rma_cel_2sd_main <- rma_cel_2sd[idx_main_2sd]
  cat("Adding HUGO genes to featuredata$geneID", fill = T)
  rma_cel_2sd_main@featureData@data$geneid <- map(rma_cel_2sd_main@featureData@data$geneassignment, split_gene_assignment_FUN) %>% sapply(unlist)
  cat("Adding ENST to featuredata$mrnaassignment", fill = T)
  rma_cel_2sd_main@featureData@data$ensgene <- map(rma_cel_2sd_main@featureData@data$mrnaassignment, split_gene_assignment_ENST_FUN) %>% sapply(unlist)
  cat("Removing unannotated probes", fill = T)
  idx_na_2sd <- which(is.na(rma_cel_2sd_main@featureData@data$geneid))
  rma_2sd <- rma_cel_2sd_main[-idx_na_2sd]
  rma_xprs <- exprs(rma_2sd)
  prb_grp_FUN <- function(goi) {
    test_probesets <- rownames(rma_2sd@featureData@data[(rma_2sd@featureData@data$geneid == goi), ])
    t_length <- length(test_probesets)
    test_xprs <- rma_xprs[which(!is.na(match(rownames(rma_xprs), test_probesets))), ]
    if (t_length == 1) {
      return(test_xprs)
    } else if (t_length == 3 | t_length == 2) {
      gen_row <- apply(test_xprs, 2, mean)
    } else if (t_length > 3) {
      gen_row <- apply(test_xprs, 2, median)
    }
  }
  genes_unq <- unique(rma_2sd@featureData@data$geneid)
  genes_unq_idx <- which(duplicated(rma_2sd@featureData@data$geneid) == FALSE)
  cat("Reducing probesets to genesets", fill = T)
  rma_xprs_grpd <- t(sapply(genes_unq, prb_grp_FUN))
  rma_2sd <- rma_2sd[genes_unq_idx]
  r_names_matched <- match(rownames(rma_xprs_grpd), rma_2sd@featureData@data$geneid)
  rownames(rma_xprs_grpd) <- rma_2sd@featureData@data$probesetid[r_names_matched]
  colnames(rma_xprs_grpd) <- colnames(rma_xprs)
  exprs(rma_2sd) <- rma_xprs_grpd
  saveRDS(rma_2sd, "./RDS/rma_2sd.rds")
  rma_2sd
}
