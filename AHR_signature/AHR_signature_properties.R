# This file explains some properties of the AHR signature genes

## Load libraries
library(purrr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(dplyr)
library(psych)

# Read the AHR signature file
overlapping_genes <- read.delim("./Signature/overlapping_AHR_signature.txt", sep = "\t", stringsAsFactors = F)

# tcga_names
tcga_names <- c("TCGA_ACC", "TCGA_BLCA", "TCGA_BRCA", "TCGA_CESC", "TCGA_COAD", "TCGA_CHOL", "TCGA_COAD", "TCGA_DLBC",
                "TCGA_ESCA", "TCGA_GBM", "TCGA_HNSC", "TCGA_KICH", "TCGA_KIRC", "TCGA_KIRP", "TCGA_LGG", "TCGA_LIHC",
                "TCGA_LUAD", "TCGA_LUSC", "TCGA_MESO", "TCGA_OV", "TCGA_PAAD", "TCGA_PCPG", "TCGA_PRAD", "TCGA_READ",
                "TCGA_SARC", "TCGA_SKCM", "TCGA_STAD", "TCGA_TGCT", "TCGA_THCA", "TCGA_THYM", "TCGA_UCEC", "TCGA_UCS", "TCGA_UVM")

# TCGA_counts (log2(cpm+1))
TCGA_counts <- readRDS("./RDS/TCGA_counts.rds")
names(TCGA_counts) <- tcga_names

# TCGA_GSVA
TCGA_gsva <- readRDS("./RDS/TCGA_GSVA_scores_safely.rds")
# Extract the GSVA scores from the safely run output
TCGA_GSVA <- TCGA_gsva[tcga_names]
TCGA_GSVA <- map(TCGA_GSVA, function(x){x$result})

###########################################################################
## Define the collinearity between the genes across the different tumors ##
###########################################################################

# Extract the counts of the AHR genes from the different tumors
AHR_cpm_genes <- map(TCGA_counts, function(x, y){x[rownames(x) %in% overlapping_genes$Gene,]})
names(AHR_cpm_genes) <- tcga_names

## Correlate between the AHR genes
AHR_signature_own_corrs <- map(AHR_cpm_genes, function(x){Hmisc::rcorr(t(x))})

# Plot the correlation of the AHR activation genes to one another ####
pdf("./Figures/AHR_all_corplots.pdf", width = 16, height = 8, pointsize = 12)
par(mfrow=c(4,8))
map2(AHR_signature_own_corrs, names(AHR_signature_own_corrs),function(cm, ct){
  corrplot(corr=cm$r,p.mat = cm$P,title = gsub("TCGA_","",ct),mar=c(0,0,1,0),cl.pos = "b",
           outline = FALSE,addgrid.col = NA, col = colorRampPalette(colors = c("blue", "white", "red"))(99), 
           bg=NULL,method = "circle",insig = "blank", order = "hclust", tl.pos = "n")
})
dev.off()

################################################################################################################
## Which AHR genes are correlating more positively or negatively to the remaining AHR signature target genes? ##
################################################################################################################

AHR_signature_own_corrs_types <- map2(AHR_signature_own_corrs, names(AHR_signature_own_corrs),function(cormat, cname){
  corrs_tables_r <- cormat$r
  corrs_tables_p <- cormat$P
  res <- map(1:163, function(i){
    df_m <- as.data.frame(cbind(corrs_tables_r[,i], corrs_tables_p[,i]))
    df_m$r <- 0
    df_m$r[which(df_m[,2]<=0.05)] <- df_m[,1][which(df_m[,2]<=0.05)]
    df_m$r
  }) %>% do.call("cbind",.) %>% as.data.frame()
  colnames(res) <- colnames(corrs_tables_p)
  rownames(res) <- colnames(corrs_tables_p)
  df <- data.frame(Tumor=cname, pos=rowSums(res>0), neg =rowSums(res<0), stringsAsFactors = F)
  df <- df %>% mutate(pos_rank = dense_rank((pos)))
  df <- df %>% mutate(neg_rank = dense_rank(desc(neg)))
  df
})

AHR_signature_own_corrs_pos_mat <- map(AHR_signature_own_corrs_types, function(x)x$pos_rank) %>% do.call(rbind,.)
colnames(AHR_signature_own_corrs_pos_mat) <- rownames(AHR_signature_own_corrs$TCGA_ACC$r)

AHR_signature_own_corrs_neg_mat <- map(AHR_signature_own_corrs_types, function(x)x$neg_rank) %>% do.call(rbind,.)
colnames(AHR_signature_own_corrs_neg_mat) <- rownames(AHR_signature_own_corrs$TCGA_ACC$r)

cor_rank_df <- data.frame(pos=apply(AHR_signature_own_corrs_pos_mat,2,sum),
                          neg= apply(AHR_signature_own_corrs_neg_mat,2,sum), stringsAsFactors = F)

cor_rank_df <- cor_rank_df[order(-cor_rank_df$pos, -cor_rank_df$neg),]
cor_rank_df$to_col <- c(rep(c("red"),10), rep("black", 153)) %>% as.factor()

# This scatter plot shows the number of positive corrs and neg corrs based ranking of the AHR genes one to another.
cor_rank_df_p <- ggscatter(cor_rank_df, x="pos", y="neg", repel = T, color = "to_col",size = 2,font.label = c(18),palette = "Dark2",rug = TRUE,legend = "none",
                           label = rownames(cor_rank_df), label.select = rownames(cor_rank_df)[1:10])
ggsave("./Figures/AHR_genes_cor_rank_df_p_plot.jpeg", device = "jpeg",plot = cor_rank_df_p, dpi = 600, height = 8, width = 12, units = "in")

#############################################################################
## Show the distribution of the GSVA scores and compare them across tumors ##
#############################################################################
## Plot GSVA score and observe their distribution
GSVA_melt_df <- map(TCGA_GSVA, t) %>% 
  map2(.,names(.), function(x,y){data.frame(Tumor=y, GSVA=x, stringsAsFactors=FALSE)}) %>%
  do.call(rbind,.)

GSVA_dist_p <- ggdensity(GSVA_melt_df, x="GSVA", facet.by = "Tumor",fill = "Tumor", color = "Tumor", add = "median", legend="none")

ggsave("./Figures/GSVA_density_plot.jpeg", device = "jpeg",plot = GSVA_dist_p, dpi = 600, height = 10, width = 10, units = "in")

############################################################################################################
## Create a correlation heatmap of the AHR signature genes with AHR activation (GSVA score) in all tumors ##
############################################################################################################

# Extract the gene signature expression values from the expression matrices
AHR_genes_gsva_dfs <- map2(TCGA_GSVA, AHR_cpm_genes, function(gsva, dge, gois){
  df <- data.frame(gsva_sc=as.numeric(gsva[1,]), t(dge[rownames(dge) %in% gois,]), stringsAsFactors = F)
  df},gois=rownames(AHR_cpm_genes$TCGA_ACC))

# pearson correlation of GSVA and AHR signature genes
AHR_genes_gsva_dfs_corrs <- map(AHR_genes_gsva_dfs, function(x){
  Hmisc::rcorr(x = as.matrix(x[,-1]), y=as.numeric(x[,1]))
})

# View the correlations of AHR activation and the AHR signature genes in a corr_heatmap
AHR_genes_gsva_dfs_corrs_tables_r <- map(AHR_genes_gsva_dfs_corrs, function(x,z){
  eois <- x$r[match(z,rownames(x$r))%>%.[!is.na(.)],"y"]
  all_eois <- vector("numeric", length = length(z))
  names(all_eois) <- z
  all_eois[match(names(eois),names(all_eois))%>%.[!is.na(.)]] <- eois
  all_eois
},z=rownames(AHR_cpm_genes$TCGA_ACC)) %>% do.call("rbind",.) %>% t()

AHR_genes_gsva_dfs_corrs_tables_p <- map(AHR_genes_gsva_dfs_corrs, function(x,z){
  eois <- x$P[match(z,rownames(x$P))%>%.[!is.na(.)],"y"]
  all_eois <- vector("numeric", length = length(z))
  names(all_eois) <- z
  all_eois[match(names(eois),names(all_eois))%>%.[!is.na(.)]] <- eois
  all_eois
},z=rownames(AHR_cpm_genes$TCGA_ACC)) %>% do.call("rbind",.) %>% t()

pdf("./Figures/corr_plot_AHRgenes_vs_GSVA_scores_AHR.pdf",width = 12, height = 8, pointsize = 12)
corrplot(corr = t(AHR_genes_gsva_dfs_corrs_tables_r), p.mat = t(AHR_genes_gsva_dfs_corrs_tables_p), insig = "blank",tl.cex = 0.5,
         tl.col = "black",cl.length = 6, cl.pos = "b",cl.cex = 0.5,
         is.corr = F,col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17),method = "circle")
dev.off()

################################################################################################################
## Plot a barplot of ranks of the AHR genes according to how strongly they correlate with the AHR GSVA scores ##
################################################################################################################

# This shoud give us an idea, which genes are better representatives in the final signature
AHR_genes_gsva_dfs_ranks <- map(1:163, function(i){
  df_m <- as.data.frame(cbind(AHR_genes_gsva_dfs_corrs_tables_r[i,], AHR_genes_gsva_dfs_corrs_tables_p[i,]))
  df_m$r <- 0
  df_m$r[which(df_m$V2<=0.05)] <- df_m$V1[which(df_m$V2<=0.05)]
  df_m$r
}) %>% do.call("rbind",.) %>% as.data.frame()

colnames(AHR_genes_gsva_dfs_ranks) <- colnames(AHR_genes_gsva_dfs_corrs_tables_r)
rownames(AHR_genes_gsva_dfs_ranks) <- rownames(AHR_genes_gsva_dfs_corrs_tables_r)

# rank based on non significant correlations
AHR_genes_gsva_dfs_ranks$zeros <- rowSums(AHR_genes_gsva_dfs_ranks == 0)

# rank based on positive significant correlations
AHR_genes_gsva_dfs_ranks$pos <- rowSums(AHR_genes_gsva_dfs_ranks > 0)-1

# rank based on negative significant correlations
AHR_genes_gsva_dfs_ranks$neg <- rowSums(AHR_genes_gsva_dfs_ranks < 0)

AHR_genes_gsva_dfs_ranks <- AHR_genes_gsva_dfs_ranks %>% mutate(z_rank = dense_rank(desc(zeros)))
AHR_genes_gsva_dfs_ranks <- AHR_genes_gsva_dfs_ranks %>% mutate(p_rank = dense_rank((pos)))
AHR_genes_gsva_dfs_ranks <- AHR_genes_gsva_dfs_ranks %>% mutate(n_rank = dense_rank(desc(neg)))

# rank based on the product of the average positive correlation and the rank of significant positive correlations
AHR_genes_gsva_dfs_ranks_fischer <- fisherz(AHR_genes_gsva_dfs_ranks[,1:32])
AHR_genes_gsva_dfs_pos_means <- apply(AHR_genes_gsva_dfs_ranks_fischer,1, function(x){mean(x[which(x>0)])}) %>% fisherz2r()
AHR_genes_gsva_dfs_ranks$r_m <- AHR_genes_gsva_dfs_pos_means*AHR_genes_gsva_dfs_ranks$p_rank
AHR_genes_gsva_dfs_ranks$r_m[is.nan(AHR_genes_gsva_dfs_ranks$r_m)] <- 0
AHR_genes_gsva_dfs_ranks <- AHR_genes_gsva_dfs_ranks %>% mutate(r_rank = dense_rank(r_m))

# Sum of ranks
AHR_genes_gsva_dfs_ranks <- AHR_genes_gsva_dfs_ranks %>% mutate(sum_of_ranks = (z_rank+p_rank+n_rank+r_rank))
AHR_genes_gsva_dfs_ranks <- AHR_genes_gsva_dfs_ranks %>% mutate(rp_rank = dense_rank(desc(sum_of_ranks)))
rownames(AHR_genes_gsva_dfs_ranks) <- rownames(AHR_genes_gsva_dfs_corrs_tables_r)

# Order the dataframe based on the sum of ranks
AHR_genes_gsva_dfs_ranks <- AHR_genes_gsva_dfs_ranks[order(-AHR_genes_gsva_dfs_ranks$sum_of_ranks),]
AHR_genes_gsva_dfs_ranks$plot <- factor(1:163)

# to avoid the overlapping of gene names
avoid_overlap <- function(x){
  ind <- seq_along(x) %% 2==0
  x[ind] <- gsub("\\s", " ", format(x[ind], width=max(nchar(x[ind]))+8))
  x[ind] <- paste0(x[ind], "                    ")
  x
}

AHR_genes_gsva_ranks_p <- ggplot(AHR_genes_gsva_dfs_ranks, aes(x=plot, y=sum_of_ranks, fill=sum_of_ranks))+geom_bar(stat = "identity",width = 0.6, position = position_dodge(width = 0.4))+
  theme(axis.text.x = element_text(aes(label=rownames(AHR_genes_gsva_dfs_ranks))))+theme_bw()+scale_fill_continuous(low = "blue", high = "red")+
  theme(axis.text.x=element_text(size = 12))+xlab("Tumors")+ylab("Sum of Ranks")+theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(breaks = 1:163, labels = avoid_overlap(rownames(AHR_genes_gsva_dfs_ranks)))

ggsave("./Figures/AHR_genes_gsva_ranks_p_plot.jpeg", device = "jpeg",plot = AHR_genes_gsva_ranks_p, dpi = 600, height = 10, width = 18, units = "in")
