#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript/AHR_signature/
## Date: Oct 2018
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This files describes the AHR signature gene ontologies
################################################

## Load libraries
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

## Read the AHR signature file (this file is already supplied for faster runs with no need for regenration)
overlapping_genes <- read.delim("../Resources/overlapping_AHR_signature_genes.txt", sep = "\t", stringsAsFactors = F)

##################################################
# Gene ontologies enriched in the AHR signature ##
##################################################
# GO analysis
GO_ALL <- map(c("BP","MF","CC"), function(x){enrichGO(gene = overlapping_genes$EID, OrgDb =org.Hs.eg.db, ont = x,
                                                      pAdjustMethod = "bonferroni",readable = TRUE, pvalueCutoff = 0.01)})
names(GO_ALL) <- c("BP","MF","CC")

# Use the GO_SemSim algorithm to group similar ontology terms 
GO_ALL_simple <- map(GO_ALL,simplify)

# Filter the higher levels
GO_ALL_filtered <- map(GO_ALL_simple,gofilter)

saveRDS(GO_ALL, "../Results/RDS/GO_all.rds")
saveRDS(GO_ALL_simple, "../Results/RDS/GO_all_simple.rds")
saveRDS(GO_ALL_filtered, "../Results/RDS/GO_all_filtered.rds")

## The result of the filtered ontologies were saved as a data.frame and were grouped manually
## Generating a plot of GO of AHR signature genes after the manual grouping of terms
go_to_plot_BP <- read.delim("../Resources/GO_table_AHR_Signatures.csv", sep = "\t", stringsAsFactors = F)
go_to_plot_BP$p.adjust <- gsub(",","\\.",go_to_plot_BP$p.adjust) %>% as.numeric()
go_to_plot_BP$logpv <- -log10(go_to_plot_BP$p.adjust)
go_to_plot_BP$GO_groups <- factor(go_to_plot_BP$GO_groups)

empty_bar <- 4
to_add <- data.frame(matrix(NA, empty_bar*nlevels(go_to_plot_BP$GO_groups), ncol(go_to_plot_BP)) )
colnames(to_add) <- colnames(go_to_plot_BP)
to_add$GO_groups <- rep(levels(go_to_plot_BP$GO_groups), each=empty_bar)
go_to_plot_BP2 <- rbind(go_to_plot_BP, to_add)
go_to_plot_BP2 <- go_to_plot_BP2 %>% arrange(GO_groups, desc(logpv))
go_to_plot_BP2$pid <- seq(1, nrow(go_to_plot_BP2))

label_data <- go_to_plot_BP2
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$pid-0.5) /number_of_bar
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

base_data <- go_to_plot_BP2 %>%
  group_by(GO_groups) %>%
  summarize(start=min(pid), end=max(pid) - empty_bar) %>%
  rowwise() %>%
  mutate(title=mean(c(start, end)))

p_GO_counts <- ggplot(go_to_plot_BP2, aes(x=as.factor(pid), y=Count, fill=logpv, group=GO_groups))+
  geom_bar(stat="identity", position = position_stack(reverse = TRUE), alpha=0.5)+ ylim(c(-30, 40))+
  theme_minimal()+ xlab("BP")+
  scale_fill_gradientn(colours = colorRampPalette(c("#1b9e77","#d95f02","#7570b3"))(99))+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+ coord_polar()+
  geom_text(data=label_data, aes(x=pid, y=36, label=Count, hjust=hjust), color="black", fontface="bold",alpha=0.6, angle= label_data$angle, inherit.aes = FALSE )+
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5, color=GO_groups), alpha=0.8, size=2.6 , inherit.aes = FALSE, show.legend = TRUE) + 
  scale_colour_brewer(palette = "Dark2", guide=guide_legend(title = "Biological processes GOs"))

pdf("../Results/AHR_signature/AHR_GOs_circular_plot.pdf", width = 8, height = 8, pointsize = 18)
p_GO_counts
dev.off()
