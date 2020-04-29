# This file describes the validation of the AHR signature using publicaly available datasets
# that were not used as part of the signature generation.

## Load libraries
library(ggplot2)
library(RColorBrewer)

## Read the AHR signature file
overlapping_genes <- read.delim("./Signature/overlapping_AHR_signature.txt", sep = "\t", stringsAsFactors = F)

####################################
## Opitz Nature 2011, Kyn_U87_8Hr ##
####################################
kyn_8h <- read.delim("./Validation/opitz_nature/8h.csv")
kyn_8h_df <- kyn_8h[which(duplicated(kyn_8h$X)==FALSE),]
kyn_8h_df <- kyn_8h_df[order(kyn_8h_df$X.2, decreasing = T),]
kyn_8h_df$AHR <- "no"
kyn_8h_df$AHR[match(overlapping_genes$Gene, kyn_8h_df$X)] <- "yes"
kyn_8h_df <- kyn_8h_df[order(kyn_8h_df$X.2,kyn_8h_df$AHR, decreasing = T),]

# plotting AHR targets only
kyn_8h_df_AHR <- kyn_8h_df[kyn_8h_df$AHR=="yes",]
kyn_8h_df_AHR$X <- factor(kyn_8h_df_AHR$X)
kyn_8h_df_AHR$ord <- factor(order(1:length(kyn_8h_df_AHR$X), decreasing = T))
labss <- kyn_8h_df_AHR$X[length(kyn_8h_df_AHR$X):1]

# plotting AHR targets in the top 10 genes
kyn_8h_df2 <- kyn_8h_df[1:10,]
kyn_8h_df2$X <- factor(kyn_8h_df2$X)
kyn_8h_df2$ord <- factor(10:1)
kyn_8h_df2$AHR <- relevel(factor(kyn_8h_df2$AHR), "yes")
labss2 <- kyn_8h_df2$X[10:1]

p_opitz_8h_t10 <- ggplot(kyn_8h_df2, aes(x=factor(kyn_8h_df2$ord), y=X.2, label=X.2, color=AHR)) + 
  geom_point(stat='identity', fill="black", size=5)+coord_flip()+
  geom_segment(aes(y = 0, x = ord, yend = X.2, xend = ord),color = "black") +ylab("\n Fold Change")+xlab("")+scale_x_discrete(labels=labss2)+scale_color_brewer(palette = "Set1")+theme_bw()+
  theme(axis.text.y = element_text(size = 14),axis.title.x = element_text(size = 14,face = "plain"), axis.text.x = element_text(size=14), legend.text = element_text(size = 14), panel.grid = element_blank())

ggsave("./Figures/Opitz_plot.jpeg", device = "jpeg",plot = p_opitz_8h_t10, dpi = 600, height = 10, width = 10, units = "in")

## Analyze additional arrays - Pending ####