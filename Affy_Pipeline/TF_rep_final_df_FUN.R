## This function defines the representation of different transcription factor target genes among differentially regulated genes
TF_rep_final_df_FUN <- function (t_table, rep_df, fc_filt, res_path, wd = 6, ht = 6, pix = 300, ...){
  matched_tt_TFs <- match(gsub("_.*", "", rownames(rep_df)), as.character(t_table$hGene))
  comp_numbers_comparison <- data.frame(TFs = rownames(rep_df), noFC = rep_df[, 1],
                                        FC = rep_df[, 2], pv = rep_df[, 3], apv = rep_df[, 4],
                                        tt_exp = t_table[matched_tt_TFs, "logFC"], group = 0)
  comp_numbers_comparison_expTFs <- comp_numbers_comparison[!is.na(comp_numbers_comparison$tt_exp), ]
  comp_idx <- which(comp_numbers_comparison_expTFs$tt_exp > 0 & comp_numbers_comparison_expTFs$tt_exp < fc_filt)
  comp_idx2 <- which(comp_numbers_comparison_expTFs$tt_exp < 0)
  comp_idx3 <- which(comp_numbers_comparison_expTFs$tt_exp >= fc_filt)
  comp_numbers_comparison_expTFs$group[comp_idx] <- "average"
  comp_numbers_comparison_expTFs$group[comp_idx2] <- "low"
  comp_numbers_comparison_expTFs$group[comp_idx3] <- "high"
  p <- ggplot(comp_numbers_comparison_expTFs, aes(x = 1:dim(comp_numbers_comparison_expTFs)[[1]], y = tt_exp))+
    geom_point(aes(col = group), stat = "identity")+ theme_bw()+ xlab("TFs")+ ylab("log2FC")
  ggsave(filename = paste0(getwd(), "/Figures/", "TFs_xprsn_distro.jpg"), plot = p, width = wd,
         height = ht, dpi = pix, units = "in")
  cat("TFs expression distribution image saved to file", fill = T)
  write.table(comp_numbers_comparison_expTFs, paste0(getwd(), res_path, "comp_numbers_comparison_TFs.txt"),
              row.names = F, quote = F, sep = "\t")
  comp_numbers_comparison_expTFs
}
