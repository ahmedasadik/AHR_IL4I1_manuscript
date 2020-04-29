## This function is a wrapper around the camera function in the limma package
gsa_FUN <- function (top_t = NULL, rma_obj, exprmnt, cntrst = "simple", 
                     msig.data.lists, cnt_mat = NULL, d_m = NULL, msigs = NULL, 
                     wts = NULL, pltfrm = "AFFY", res_path){
  gsa_fun <- function(cont, idcs, Mmat, des, weit, msig_l){
    gsa <- map2(idcs, rep(list(cont), by = msig_l), function(a, b) {
      camera(y = Mmat, index = a, design = des, contrast = b, weights = weit, trend.var = T)})
    names(gsa) <- paste(names(cont)[which(cont == 1)], names(idcs), sep = "_")
    gsa
  }
  platform <- tolower(pltfrm)
  if (platform == "affy") {
    cat("Affy Hugene-2.0-st platform selected", fill = TRUE)
    prbsets_gsa_2sd <- rownames(top_t)
    gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_obj@featureData@data$probesetid)
    gsa_2sd_eset <- rma_obj[gsa_2sd_idx]
    gsa_all_mapped_idcs <- map(msig.data.lists, ids2indices, top_t$hEID)
    cat("Mapping idcs completed", fill = TRUE)
    if (cntrst == "simple") {
      cat("Simple comparison ctrl_trt selected", fill = TRUE)
      cat("Running Camera", fill = TRUE)
      gsa_all_res_2sd <- map(gsa_all_mapped_idcs, function(id) {
        camera(y = exprs(gsa_2sd_eset), index = id, design = d_m, contrast = 2, trend.var = T)})
      walk2(gsa_all_res_2sd, paste(res_path, "gsa_", exprmnt, "_res_", names(gsa_all_mapped_idcs), ".txt", sep = ""), write.table, sep = "\t")
      cat("GSA comparison files written to disk", fill = TRUE)
      return(gsa_all_res_2sd)
    } else {
      cont_names <- gsub(" \\- ", "_", colnames(cnt_mat))
      if (length(cont_names) > 1) {
        cat("Multiple contrasts selected", fill = TRUE)
        cnt_list <- list()
        for (i in seq_along(cont_names)) {
          cnt_list[[i]] <- cnt_mat[, i]
        }
        names(cnt_list) <- cont_names
        cat("List of contrasts created", fill = TRUE)
        no_cores <- detectCores() - 2
        cl <- makeCluster(no_cores)
        clusterEvalQ(cl, {library(limma);library(purrr)})
        cat("Running Camera", fill = TRUE)
        gsa_cont <- parLapplyLB(cl = cl, X = cnt_list, 
                                fun = gsa_fun, idcs = gsa_all_mapped_idcs, 
                                Mmat = exprs(gsa_2sd_eset), des = d_m, msig_l = msigs, 
                                weit = wts)
        stopCluster(cl)
        for (i in seq_along(cont_names)) {
          walk2(gsa_cont[[i]], paste(res_path, "gsa_", cont_names[i], "_", exprmnt, "_res_", names(gsa_all_mapped_idcs), ".txt", sep = ""), write.table, sep = "\t")
        }
        cat("GSA comparison files written to disk", fill = TRUE)
        return(gsa_cont)
      } else {
        gsa_all_res_2sd <- map(gsa_all_mapped_idcs, function(id) {
          camera(y = exprs(gsa_2sd_eset), index = id, design = d_m, contrast = cnt_mat[, 1], trend.var = T)})
        walk2(gsa_all_res_2sd, paste(res_path, "gsa_", exprmnt, "_res_", names(gsa_all_mapped_idcs), ".txt", sep = ""), write.table, sep = "\t")
        cat("GSA comparison files written to disk", fill = TRUE)
        return(gsa_all_res_2sd)
      }
    }
  }
  else {
    cat(paste0("Platform ", platform, " is selected"), fill = TRUE)
    gsa_all_mapped_idcs <- map(msig.data.lists, ids2indices, rma_obj$genes$EID)
    cat("Mapping idcs completed", fill = TRUE)
    if (cntrst == "simple") {
      cat("Simple comparison ctrl_trt selected", fill = TRUE)
      cat("Running Camera", fill = TRUE)
      gsa_all_res_2sd <- map(gsa_all_mapped_idcs, function(id) {
        camera(y = rma_obj$M, index = id, design = d_m, contrast = 2, trend.var = T)})
      walk2(gsa_all_res_2sd, paste(res_path, "gsa_", exprmnt, "_res_", names(gsa_all_mapped_idcs), ".txt", sep = ""), write.table, sep = "\t")
      cat("GSA comparison files written to disk", fill = TRUE)
      return(gsa_all_res_2sd)
    } else {
      cont_names <- gsub(" \\- ", "_", colnames(cnt_mat))
      if (length(cont_names) > 1) {
        cat("Multiple contrasts selected", fill = TRUE)
        cnt_list <- list()
        for (i in seq_along(cont_names)) {
          cnt_list[[i]] <- cnt_mat[, i]
        }
        names(cnt_list) <- cont_names
        cat("List of contrasts created", fill = TRUE)
        no_cores <- detectCores() - 2
        cl <- makeCluster(no_cores)
        clusterEvalQ(cl, {library(limma);library(purrr)})
        cat("Running Camera", fill = TRUE)
        gsa_cont <- parLapplyLB(cl = cl, X = cnt_list, fun = gsa_fun, idcs = gsa_all_mapped_idcs,
                                Mmat = rma_obj$M, des = d_m, msig_l = msigs, weit = wts)
        stopCluster(cl)
        tt <- map2(gsa_cont, cont_names, function(a, b, xprmnt, go_names, res_path = res_path) {
          walk2(a, paste(res_path, "gsa_", b, "_", xprmnt, "_res_", go_names, ".txt", sep = ""), write.table, sep = "\t")
        }, xprmnt = exprmnt, go_names = names(gsa_all_mapped_idcs))
        cat("GSA comparison files written to disk", fill = TRUE)
        return(gsa_cont)
      } else {
        gsa_all_res_2sd <- map(gsa_all_mapped_idcs, function(id) {
          camera(y = rma_obj$M, index = id, design = d_m, contrast = cnt_mat[, 1], trend.var = T)})
        walk2(gsa_all_res_2sd, paste(res_path, "gsa_", exprmnt, "_res_", names(gsa_all_mapped_idcs), ".txt", sep = ""), write.table, sep = "\t")
        cat("GSA comparison files written to disk", fill = TRUE)
        return(gsa_all_res_2sd)
      }
    }
  }
}