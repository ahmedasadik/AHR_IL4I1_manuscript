## This function is a wrapper of the camera function used for the results of transcription factor analysis
TF_camera_gsa_FUN <- function(TF_list, top_T = NULL, rma_obj, pltfrm = "Affy", d_m, 
                              cntrst = NULL, cnt_mat = NULL, TF_motifs = FALSE, wts = NULL, 
                              res_path){
  if (pltfrm == "Affy") {
    print("Affy hugene-2.0-st selected")
    prbsets_gsa_2sd <- rownames(top_T)
    gsa_2sd_idx <- match(prbsets_gsa_2sd, rma_obj@featureData@data$probesetid)
    gsa_2sd_eset <- rma_obj[gsa_2sd_idx]
    TFSBDB_list_idcs <- ids2indices(TF_list, top_T$hEID)
    TFSBDB_list_enrichment <- camera(y = exprs(gsa_2sd_eset), index = TFSBDB_list_idcs,
                                     design = d_m, contrast = cntrst, weights = wts, trend.var = T)
  } else {
    print("Other platform selected")
    TFSBDB_list_idcs <- ids2indices(TF_list, rma_obj$genes$EID)
    TFSBDB_list_enrichment <- camera(y = rma_obj, index = TFSBDB_list_idcs, design = d_m,
                                     contrast = cntrst, weights = wts, trend.var = T)
  }
  if (TF_motifs == TRUE) {
    write.table(TFSBDB_list_enrichment, paste(res_path, "TFSBDB_list_enrichment_separate_motifs.txt", sep = ""), sep = "\t")
  } else {
    write.table(TFSBDB_list_enrichment, paste(res_path, "TFSBDB_list_enrichment_combined_motifs.txt", sep = ""), sep = "\t")
  }
  return(TFSBDB_list_enrichment)
}
