#!/usr/bin/env Rscript

#################################################
## Project: AHR mouse microarray analysis
## Date: Oct 2018
## Author: Murat Iskar
##
## Description
## This script describes how we analyze mouse microarray samples of effector and memory cells from IL4I1 knock-out and wild type mouse.
#################################################


options(max.print = 200)
options(stringsAsFactors = FALSE)

library(affy)
library(ggplot2)
library(affycoretools)
library(ggsci)
library(corrplot)
library(limma)
library(biomaRt)
library(gplots)
library(mogene20sttranscriptcluster.db)
library(heatmap3)

################################################################################
################################################################################
################################################################################
#Normalization of all the microarray samples.
cels = list.files(
  "../raw_data/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/affy_chips01",
  pattern = "CEL",
  full.names = TRUE
)
all.eset = oligo::rma(oligo::read.celfiles(cels))
all.eset <-
  annotateEset(
    all.eset,
    mogene20sttranscriptcluster.db,
    columns(mogene20sttranscriptcluster.db)
  )
m = exprs(all.eset)
k = featureData(all.eset)

write.table(m, file = "../results/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-rma-matrix.txt", 
				sep = "\t")

all.eset = getMainProbes(all.eset, level = "core")

################################################################################
#Filter lowly expressed genes:
m = exprs(all.eset)
all_medians <- rowMedians(Biobase::exprs(all.eset))

#Checking the signal distribution of probe sets.
pdf(
  paste0(
    "../results/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/figures/",
    Sys.Date(),
    "-MI-median-plot.pdf"
  )
)
hist(
  all_medians,
  100,
  col = "black",
  freq = FALSE,
  main = "Histogram of the median intensities",
  border = "antiquewhite4",
  xlab = "Median intensities"
)
dev.off()
#threshold set to 4.
#median of any category should be higher than 4.
select = (rowMedians(m[, grep("Effector KO", colnames(m))]) > 4 |
            rowMedians(m[, grep("Effector WT", colnames(m))]) > 4 |
            rowMedians(m[, grep("Memory KO", colnames(m))]) > 4 |
            rowMedians(m[, grep("Memory WT", colnames(m))]) > 4)

all.eset2 <- subset(all.eset, select)

#Remove probes with missing gene symbol:
all.eset2 <-
  subset(all.eset2,!is.na(featureData(all.eset2)@data$SYMBOL))

#select probeset with highest variance per gene:
allsd = apply(Biobase::exprs(all.eset2), 1, sd)
sy = featureData(all.eset2)@data$SYMBOL
sy = sy[order(allsd, decreasing = TRUE)]
names(allsd) = featureData(all.eset2)@data$PROBEID
allsd = sort(allsd, decreasing = TRUE)
allsd = allsd[!duplicated(sy)]

all.eset2 <-
  subset(all.eset2, (featureData(all.eset2)@data$PROBEID %in% names(allsd)))
################################################################################

###############################################################################
#Principal component analysis and the visualization
expmat = exprs(all.eset2)
expmat2 = expmat[apply(expmat, 1, sd) > sort(apply(expmat, 1, sd), 
					decreasing = TRUE)[2001], ]

PCA_raw <- prcomp(t(expmat2), scale = FALSE)
prop = PCA_raw$sdev ^ 2 / sum(PCA_raw$sdev ^ 2)
dataGG <- data.frame(PC1 = PCA_raw$x[, 1], PC2 = PCA_raw$x[, 2])
colsidecol = colnames(expmat2)
colsidecol[grep("Effector KO", colsidecol)] = "eff KO"
colsidecol[grep("Effector WT", colsidecol)] = "eff WT"
colsidecol[grep("Memory KO", colsidecol)] = "me KO"
colsidecol[grep("Memory WT", colsidecol)] = "me WT"
Ydf = data.frame(x = PCA_raw$x[, 1],
                 y = PCA_raw$x[, 2],
                 colsidecol = colsidecol)


#Preparing PCA figure in Supplementary Figure 7D
ggplot(Ydf, aes(x = x, y = y, col = colsidecol)) + geom_point(
  aes(fill = colsidecol),
  size = 17,
  colour = "black",
  pch = 21
) + theme_bw() +
  scale_color_npg() + scale_fill_npg() + theme_bw() + theme(panel.spacing = unit(0, "lines")) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  xlab(paste0("PC1 (", round(prop[1] * 100, 0), "%)")) + ylab(paste0("PC2 (", round(prop[2] *
                                                                                      100, 0), "%)")) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

ggsave(
  file = paste0(
    "../results/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/figures/",
    Sys.Date(),
    "-MI-PCA-plot-all-genesfiltered-noNA-lowexpfiltered-Supp-fig-7D.pdf"
  ),
  width = 10,
  height = 7,
  useDingbats = FALSE
)

###############################################################################
#Plotting Heatmap:
#Preparing Figure 6F and Supplementary Figure 7E

expmat = exprs(all.eset2)

k = abs(apply(expmat[, 1:4], 1, mean) - apply(expmat[, 5:8], 1, mean))
effsel = which(k > sort(k, decreasing = TRUE)[501])

combexp = expmat[effsel, ]
#Median centering
combexp = cbind(combexp[, 1:8] - apply(combexp[, 1:8], 1, median),
                combexp[, 9:16] - apply(combexp[, 9:16], 1, median))

colsidecol = colnames(combexp)
colsidecol[grep("Effector KO", colsidecol)] = "eff KO"
colsidecol[grep("Effector WT", colsidecol)] = "eff WT"
colsidecol[grep("Memory KO", colsidecol)] = "mem KO"
colsidecol[grep("Memory WT", colsidecol)] = "mem WT"
colnames(combexp) = colsidecol

combexp = combexp[!is.na(combexp[, 1]), ]
#Genes to be highlighted:
highgenes = t(read.table("highlight-genes-mogene-st2.txt"))


symdata = featureData(all.eset)@data

nam = rownames(combexp)
nam = symdata[match(t(nam), symdata[, 1]), ]$SYMBOL
genes = as.character(symdata[match(t(highgenes), symdata[, 1]), ]$SYMBOL)
#Adding two more genes to be highlighted
genes = c(genes, "Gzma" , "Ahr")
nam = as.character(nam)
names(nam) = nam
nam[!(nam %in% genes)] = ""

colsidecol = rep(c("green", "blue", "orange", "red"), each = 4)


h = hclust(as.dist(1 - cor(t(combexp))), method = "average")
combexp2 = combexp[h$order, ]
nam2 = nam[h$order]
combexp2[combexp2 > 1.5] = 1.5
combexp2[combexp2 < (-1.5)] = (-1.5)

addcol = as.numeric(nam2 != "")
addcol[addcol == 0] = NA
addcol[addcol == 1] = 0
combexp2 = cbind(combexp2, addcol)
colsidecol = c(colsidecol, "black")

#color ramp based on viridis package:
load("colramp.RData")


#Heatmap in pdf:
library(heatmap3)

#Ordering the samples based on the mean absolute expression change after median-centering.
#the ordering of memory samples is based on matched effector samples.
combexp2 = combexp2[, c(
  order(apply(abs(combexp2[, 1:4]), 2, sum), decreasing = TRUE),
  order(apply(abs(combexp2[, 5:8]), 2, sum), decreasing = TRUE) + 4,
  order(apply(abs(combexp2[, 1:4]), 2, sum), decreasing = TRUE) + 8,
  order(apply(abs(combexp2[, 5:8]), 2, sum), decreasing = TRUE) + 12,
  17,
  18
)]


pdf(
  paste0(
    "../results/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/figures/",
    Sys.Date(),
    "-MI-heatmap-diff-genes-in-eff-mem-ko-vs-wt-Figure6F-SuppFigure7E.pdf"
  )
)
heatmap3(
  combexp2,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  balanceColor = TRUE,
  labRow = nam2,
  col = colramp,
  ColSideColors = colsidecol
)
dev.off()
###############################################################################
#Differential analysis:

expmat = exprs(all.eset2)

treatment = rep(0, 16)
treatment[grep("Effector KO", colnames(expmat))] = 1
treatment[grep("Effector WT", colnames(expmat))] = 2
treatment[grep("Memory KO", colnames(expmat))] = 3
treatment[grep("Memory WT", colnames(expmat))] = 4
treatment <- factor(c(treatment))
design <- model.matrix( ~ 0 + treatment)

fit <- lmFit(expmat, design)

contrast.matrix <-
  makeContrasts(treatment1 - treatment2, treatment3 - treatment4, levels =
                  design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

result1 <-
  topTable(fit2,
           number = length(expmat[, 1]),
           coef = 1,
           adjust = "fdr")
result1 <- result1[order(rownames(result1)), ]

result2 <-
  topTable(fit2,
           number = length(expmat[, 1]),
           coef = 2,
           adjust = "fdr")
result2 <- result2[order(rownames(result2)), ]

allres = cbind(featureData(all.eset2)@data, result1, result2)



#Preparing supplementary table for mouse microarray differential analysis
write.table(
  allres,
  file = paste0(
    "../results/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/figures/",
    Sys.Date(),
    "-all-differential-results-Effector-Memory-SuppTable-mouse-differential-analysis.txt"
  ),
  sep = "\t",
  quote = FALSE
)


#Get orthologs mouse to human:
musGenes <- featureData(all.eset2)@data$SYMBOL

convertMouseGeneList <- function(x) {
  library("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(
    attributes = c("mgi_symbol"),
    filters = "mgi_symbol",
    values = x ,
    mart = mouse,
    attributesL = c("hgnc_symbol"),
    martL = human,
    uniqueRows = T
  )
  humanx <- unique(genesV2[, 2])
  return(genesV2)
}

humanSYMBOL = convertMouseGeneList(musGenes)
hs = humanSYMBOL[match(featureData(all.eset2)@data$SYMBOL, humanSYMBOL[, 1]), 2]
humangene = hs
allres = cbind(featureData(all.eset2)@data, humangene, result1, result2)


#Preparing supplementary table for mouse microarray differential analysis including Human gene ortohologs
write.table(
  allres,
  file = paste0(
    "../results/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/figures/",
    Sys.Date(),
    "-all-differential-results-with-human-orth.txt"
  ),
  sep = "\t",
  quote = FALSE
)

#Preparing fold change ranked list of genes (in human orthologs) to provide as input to gene-set enrichment analysis pre-ranked.
write.table(
  allres[, c(26, 27)],
  file = paste0(
    "../results/CD8effectorCD127-CD44+andmemoryCD127+CD44+MoGene-2_0-st/figures/",
    Sys.Date(),
    "-effector-FC-human-FCranked-for-GSEA.rnk"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)