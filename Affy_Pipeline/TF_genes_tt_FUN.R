## This function returns the number of TF targets present in a toptable
TF_genes_tt_FUN <- function (ptrns, t_table, fc_filt = NA, rtrn_value = "dim", pv = NA, apv = NA){
  tt_idx <- match(ptrns, t_table[, 1])
  tt_idx <- tt_idx[!is.na(tt_idx)]
  tt <- t_table[tt_idx, ]
  if (rtrn_value == "dim" && is.na(fc_filt) && is.na(pv) && is.na(apv)) {
    return(dim(tt)[[1]])
  } else if (rtrn_value == "dim" && !is.na(fc_filt) && is.na(pv) && is.na(apv)) {
    tt_fc <- tt[which(tt$logFC > fc_filt | tt$logFC < -fc_filt), ]
    return(dim(tt_fc)[[1]])
  } else if (rtrn_value == "dim" && !is.na(fc_filt) && !is.na(pv) && is.na(apv)) {
    tt_fc <- tt[which(tt$logFC > fc_filt | tt$logFC < -fc_filt), ]
    tt_pv <- tt_fc[which(tt_fc$P.Value <= pv), ]
    return(dim(tt_pv)[[1]])
  } else if (rtrn_value == "dim" && !is.na(fc_filt) && !is.na(pv) && !is.na(apv)) {
    tt_fc <- tt[which(tt$logFC > fc_filt | tt$logFC < -fc_filt), ]
    tt_pv <- tt_fc[which(tt_fc$P.Value <= pv), ]
    tt_apv <- tt_pv[tt_pv$adj.P.Val <= apv, ]
    return(dim(tt_apv)[[1]])
  } else if (rtrn_value == "tt_fc" && is.na(fc_filt) && is.na(pv) && is.na(apv)) {
    return((tt))
  } else if (rtrn_value == "tt_fc" && !is.na(fc_filt) && is.na(pv) && is.na(apv)) {
    tt_fc <- tt[which(tt$logFC > fc_filt | tt$logFC < -fc_filt), ]
    return(tt_fc)
  } else if (rtrn_value == "tt_fc" && !is.na(fc_filt) && !is.na(pv) && is.na(apv)) {
    tt_fc <- tt[which(tt$logFC > fc_filt | tt$logFC < -fc_filt), ]
    tt_pv <- tt_fc[which(tt_fc$P.Value <= pv), ]
    return(tt_pv)
  } else if (rtrn_value == "tt_fc" && !is.na(fc_filt) && !is.na(pv) && !is.na(apv)) {
    tt_fc <- tt[which(tt$logFC > fc_filt | tt$logFC < -fc_filt), ]
    tt_pv <- tt_fc[which(tt_fc$P.Value <= pv), ]
    tt_apv <- tt_pv[tt_pv$adj.P.Val <= apv, ]
    return(tt_apv)
  }
}
