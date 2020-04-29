## This function generates a barplot of the gene set testing results produced by camera
gsa_bar_plot_FUN <- function (idx, gsa_ls, ts = 8, res_path, coi, wd, ht, pix = 600){
  pthwy_name <- names(gsa_ls)[idx]
  gsa_df <- gsa_ls[[idx]]
  gsa_df$pathways <- rownames(gsa_df)
  gsa_df <- gsa_df[gsa_df$PValue <= 0.05, ]
  gsa_df <- gsa_df[gsa_df$NGenes >= 10, ]
  gsa_df$Direction <- as.factor(gsa_df$Direction)
  up_idx <- which(gsa_df$Direction == "Up")
  dn_idx <- which(gsa_df$Direction == "Down")
  up_df <- gsa_df[up_idx, ] %>% .[order(.$PValue, decreasing = F),]
  dn_df <- gsa_df[dn_idx, ] %>% .[order(.$PValue, decreasing = F),]
  if (length(up_idx) >= 10 & length(dn_idx) >= 10){
    cmb_df <- rbind(up_df[1:10, ], dn_df[1:10,])
    cmb_df$PValue <- -log10(cmb_df$PValue)
  } else if (length(up_idx) >= 10 & length(dn_idx) < 10) {
    cmb_df <- rbind(up_df[1:(10 + (10 - length(dn_idx))),], dn_df[1:length(dn_idx), ])
    cmb_df$PValue <- -log10(cmb_df$PValue)
  } else if (length(up_idx) < 10 & length(dn_idx) >= 10) {
    cmb_df <- rbind(dn_df[1:(10 + (10 - length(up_idx))),], up_df[1:length(up_idx), ])
    cmb_df$PValue <- -log10(cmb_df$PValue)
  } else {
    cmb_df <- rbind(up_df, dn_df)
    cmb_df$PValue <- -log10(cmb_df$PValue)
  }
  gsa_p <- ggplot(cmb_df, aes(x = reorder(pathways, PValue), y = PValue, fill = Direction))+
    geom_bar(stat = "identity")+ theme_bw()+
    scale_fill_manual(values = brewer.pal(8,"Set1")[c(3, 5)])+ guides(fill = guide_legend(reverse = TRUE))+ 
    ylab("-log10(PValue)") + xlab("Pathways")+
    theme(text = element_text(size = ts), axis.text.y = element_blank(),
          panel.grid = element_blank(), axis.ticks.y = element_blank())+
    ggtitle(paste0("Top ", dim(cmb_df)[[1]], " enriched ", pthwy_name, " pathways"))+ 
    geom_text(aes(label = pathways, y = 1, hjust = 0),
              position = position_nudge(0, -0.9))+
    geom_text(aes(label = NGenes, y = 1, hjust = 0), position = position_nudge(0, -1.3))+
    coord_flip()+ scale_colour_manual(values = "grey4")
  cat(toupper(pthwy_name), "saved to disk", fill = T)
  ggsave(filename = paste0(getwd(), res_path, coi, "_gsa_", pthwy_name, ".jpg"),
         plot = gsa_p, width = wd, height = ht, dpi = pix, units = "in")
}
