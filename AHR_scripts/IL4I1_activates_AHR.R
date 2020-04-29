# This script describes how we came across IL4I1 as an enzyme producing ligands activating AHR,
# in addition we also show that GBM would be the first tumor of choice for studying IL4I1

## Load libraries
library(purrr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(corrplot)
library(psych)

# tcga_names
tcga_names <- c("TCGA_ACC", "TCGA_BLCA", "TCGA_BRCA", "TCGA_CESC", "TCGA_COAD", "TCGA_CHOL", "TCGA_COAD", "TCGA_DLBC",
                "TCGA_ESCA", "TCGA_GBM", "TCGA_HNSC", "TCGA_KICH", "TCGA_KIRC", "TCGA_KIRP", "TCGA_LGG", "TCGA_LIHC",
                "TCGA_LUAD", "TCGA_LUSC", "TCGA_MESO", "TCGA_OV", "TCGA_PAAD", "TCGA_PCPG", "TCGA_PRAD", "TCGA_READ",
                "TCGA_SARC", "TCGA_SKCM", "TCGA_STAD", "TCGA_TGCT", "TCGA_THCA", "TCGA_THYM", "TCGA_UCEC", "TCGA_UCS", "TCGA_UVM")

# Selected enzymes of the Tryptophan degradation pathway
trp_enz_sel <- c("IL4I1","IDO1","TDO2","KMO","KYNU","MAOB","MAOA","CAT","TPH1","TPH2","DDC","KYAT3","KYAT1","AADAT","AFMID")

# TCGA_voom
TCGA_voom <- readRDS("./RDS/TCGA_DGE_voom_annots.rds")

# TCGA tumors and the modules correlating with AHR (this shows both positive and negative correlations - pearson)
TCGA_MEs_GSVA_cors <- readRDS("./RDS/TCGA_MEs_GSVA_cors.rds")

# Modules correlating only positively with AHR activation
modules_with_enzymes_pos <- readRDS("./RDS/modules_with_enzymes_pos.rds")

# TCGA_GSVA
TCGA_gsva <- readRDS("./RDS/TCGA_GSVA_scores_safely.rds")
# Extract the GSVA scores from the safely run output
TCGA_GSVA <- TCGA_gsva[tcga_names]
TCGA_GSVA <- map(TCGA_GSVA, function(x){x$result})

#######################################################################
## Since IDO1 and TDO2 are not driving AHR activation in all tumors, ##
## what is the status of other enzymes of Trp degradation pathway?   ##
#######################################################################

## The incidence plots of selected enzymes from Trp degradation pathway
final_TCGA_MEs_GSVA_corrs_AHR_enzymes <- map2(TCGA_MEs_GSVA_cors, modules_with_enzymes_pos, function(cors, enzs){
  cors$enz_aaas_l <- ""
  cors$enz_ahr_l <- ""
  cors$enz_aaas_g <- ""
  cors$enz_ahr_g <- ""
  cors$enz_aaas_ahr <- ""
  enz_modules <- enzs$module
  for (i in 1:length(cors$modules)){
    enz_mod_idx <- which(enz_modules==cors$modules[i])
    aaas_idx <- which(enzs[enz_mod_idx,"aaa_enzymes"]=="yes")
    ahr_idx <- which(enzs[enz_mod_idx,"AHR_genes"]=="yes")
    cors$enz_aaas_l[i] <- length(aaas_idx)
    cors$enz_ahr_l[i] <- length(ahr_idx)
    cors$enz_aaas_g[i] <- paste(enzs$genes[enz_mod_idx][aaas_idx],collapse = ";")
    cors$enz_ahr_g[i] <- paste(enzs$genes[enz_mod_idx][ahr_idx],collapse = ";")
    cors$enz_aaas_ahr[i] <- paste(enzs$genes[enz_mod_idx][intersect(aaas_idx, ahr_idx)],collapse = ";")
  }
  cors
})

names(final_TCGA_MEs_GSVA_corrs_AHR_enzymes) <- tcga_names

final_TCGA_MEs_GSVA_corrs_AHR_enzymes <- map2(final_TCGA_MEs_GSVA_corrs_AHR_enzymes, tcga_names, function(x,y){
  df <- data.frame(tumor= gsub("TCGA_","",y),x)
}) %>% do.call(rbind,.)

final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IL4I1 <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IDO1 <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$TDO2 <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$none <- ""
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IL4I1[grep("IL4I1",final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_g)] <- "yes"
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$IDO1[grep("IDO1",final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_g)] <- "yes"
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$TDO2[grep("TDO2",final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_g)] <- "yes"
final_TCGA_MEs_GSVA_corrs_AHR_enzymes$none[which(final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_l==0 & final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_ahr_l==0)] <- "ns"

enz_no_in_modules <- sapply(trp_enz_sel, function(x,y){
  length(grep(x,y))
}, y=final_TCGA_MEs_GSVA_corrs_AHR_enzymes$enz_aaas_g)

enz_no_in_modules <- data.frame(genes=names(enz_no_in_modules), hits=enz_no_in_modules)
enz_no_in_modules <- enz_no_in_modules[order(enz_no_in_modules$hits, decreasing = F),]
enz_no_in_modules$numbers <- 1:length(enz_no_in_modules$genes)
enz_no_in_modules$numbers <- as.factor(enz_no_in_modules$numbers)

avoid_overlap <- function(x){
  ind <- seq_along(x) %% 2==0
  x[ind] <- gsub("\\s", " ", format(x[ind], width=max(nchar(x[ind]))+8))
  x[ind] <- paste0(x[ind], "         ")
  x
}
a <- ifelse(as.character(enz_no_in_modules$genes) == "IL4I1", "red2","black")

tiff("./Figures/TRP_enz_ME_incidences_pos_cor.tiff",width = 5,height = 3,units = "in",res = 300)
ggplot(enz_no_in_modules, aes(x=numbers, y=hits, fill=hits))+geom_bar(stat = "identity",width = 0.6, position = position_dodge(width = 0.4))+
  theme(axis.text.x = element_text(aes(label=enz_no_in_modules$genes)))+theme_bw()+scale_fill_continuous(low = "blue", high = "red")+
  theme(axis.text.y=element_text(colour = a, size = 12))+ylab("number of modules")+xlab(label = NULL)+
  scale_x_discrete(breaks = 1:length(trp_enz_sel), labels = avoid_overlap(as.character(enz_no_in_modules$genes)))+coord_flip()
dev.off()

###############################################################################################################
## Ranking the tumors based on AHR activation correlations with the Trp degradation pathway selected enzymes ##
###############################################################################################################

## Correlations of trp enzymes and incidence rates in the modules that are positively correlating with AHR
Enz_gsva_dfs <- map2(TCGA_GSVA, TCGA_voom, function(gsva, dge, gois){
  df <- data.frame(gsva_sc=as.numeric(gsva[1,]), dge$E[,colnames(dge$E) %in% gois], stringsAsFactors = F)
  df},gois=as.character(trp_enz_sel))

# pearson correlation of GSVA and enzymes
corrs <- map(Enz_gsva_dfs, function(x){
  Hmisc::rcorr(x = as.matrix(x[,-1]), y=as.numeric(x[,1]))
})

corrs_tables_r <- map(corrs, function(x,z){
  eois <- x$r[match(z,rownames(x$r))%>%.[!is.na(.)],"y"]
  all_eois <- vector("numeric", length = length(z))
  names(all_eois) <- z
  all_eois[match(names(eois),names(all_eois))%>%.[!is.na(.)]] <- eois
  all_eois
},z=trp_enz_sel) %>% do.call("rbind",.)

corrs_tables_p <- map(corrs, function(x,z){
  eois <- x$P[match(z,rownames(x$P))%>%.[!is.na(.)],"y"]
  all_eois <- vector("numeric", length = length(z))
  names(all_eois) <- z
  all_eois[match(names(eois),names(all_eois))%>%.[!is.na(.)]] <- eois
  all_eois
},z=trp_enz_sel) %>% do.call("rbind",.)


tiff("./Figures/TRP_enz_AHR_corrs.tiff",width = 6,height = 8, units = "in", res = 300)
corrplot(corr = corrs_tables_r, p.mat = corrs_tables_p, insig = "blank",tl.cex = 0.7,
         tl.col = "black",cl.length = 6, cl.pos = "b",cl.cex = 0.5,
         is.corr = F,col = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(17),method = "circle")
dev.off()

# Creat a table of ranks for the previously made corr_heatmap of the selected Trp enzymes
Trp_enz_ranks <- map(1:32, function(i){
  df_m <- as.data.frame(cbind(corrs_tables_r[i,], corrs_tables_p[i,]))
  df_m$r <- 0
  df_m$r[which(df_m$V2<=0.05)] <- df_m$V1[which(df_m$V2<=0.05)]
  df_m$r
}) %>% do.call("rbind",.) %>% as.data.frame()

colnames(Trp_enz_ranks) <- colnames(corrs_tables_r)
rownames(Trp_enz_ranks) <- rownames(corrs_tables_r)

# rank based on non significant correlations
Trp_enz_ranks$zeros <- rowSums(Trp_enz_ranks == 0)

# rank based on positive significant correlations
Trp_enz_ranks$pos <- rowSums(Trp_enz_ranks > 0)-1

# rank based on negative significant correlations
Trp_enz_ranks$neg <- rowSums(Trp_enz_ranks < 0)

Trp_enz_ranks <- Trp_enz_ranks %>% mutate(z_rank = dense_rank(desc(zeros)))
Trp_enz_ranks <- Trp_enz_ranks %>% mutate(p_rank = dense_rank((pos)))
Trp_enz_ranks <- Trp_enz_ranks %>% mutate(n_rank = dense_rank(desc(neg)))
Trp_enz_ranks_fischer <- fisherz(Trp_enz_ranks[,1:length(trp_enz_sel)])

# rank based on the product of the average positive correlation and the rank of significant positive correlations
pos_means <- apply(Trp_enz_ranks_fischer,1, function(x){mean(x[which(x>0)])}) %>% fisherz2r()
Trp_enz_ranks$r_m <- pos_means*Trp_enz_ranks$p_rank
Trp_enz_ranks <- Trp_enz_ranks %>% mutate(r_rank = dense_rank(r_m))

# Sum of ranks
Trp_enz_ranks <- Trp_enz_ranks %>% mutate(sum_of_ranks = (z_rank+p_rank+n_rank+r_rank))
Trp_enz_ranks <- Trp_enz_ranks %>% mutate(rp_rank = dense_rank(desc(sum_of_ranks)))

# Order the dataframe based on the sum of ranks
rownames(Trp_enz_ranks) <- rownames(corrs_tables_r)
Trp_enz_ranks <- Trp_enz_ranks[order(Trp_enz_ranks$sum_of_ranks),]
Trp_enz_ranks$plot <- factor(1:32)

# to avoid the overlapping of gene names
avoid_overlap <- function(x){
  ind <- seq_along(x) %% 2==0
  x[ind] <- gsub("\\s", " ", format(x[ind], width=max(nchar(x[ind]))))
  x[ind] <- paste0(x[ind], "              ")
  x
}

tiff("./Figures/TRP_enz_tumor_ranks.tiff",width = 6,height = 4,units = "in",res = 300)
ggplot(Trp_enz_ranks, aes(x=plot, y=sum_of_ranks, fill=sum_of_ranks))+geom_bar(stat = "identity",width = 0.6, position = position_dodge(width = 0.4))+
  theme(axis.text.x = element_text(aes(label=rownames(Trp_enz_ranks))))+theme_bw()+scale_fill_continuous(low = "blue", high = "red")+
  theme(axis.text.y=element_text(size = 12))+xlab("Tumors")+ylab("Sum of Ranks")+
  scale_x_discrete(breaks = 1:32, labels = avoid_overlap(as.character(gsub("TCGA_","",rownames(Trp_enz_ranks)))))+coord_flip()
dev.off()

################################################################################################
## Scatter plot showing which tumors should be considered to study IL4I1 AHR-mediated effects ## 
################################################################################################

## Sum of Ranks vs IL4I1 corrs
tiff("./Figures/IL4I1_scatter_vs_Ranks.tiff",width = 6,height = 4,units = "in", res = 300)
ggscatter(Trp_enz_ranks, x="sum_of_ranks", y="IL4I1",color = "plot",
          repel = T, size = 3,ylab = "IL4I1 R2",
          label = gsub("TCGA_","",rownames(Trp_enz_ranks)), legend="none", font.label = c(12, "black"))
dev.off()
