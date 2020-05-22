#!/usr/bin/env Rscript

#################################################
## Project: AHR IL4I1
## Origin: https://github.com/ahmedasadik/AHR_IL4I1_manuscript
## Date: May 2020
## Author: Ahmed Sadik (a.sadik@dkfz.de)
##
## Description
## This script describes how we came across IL4I1 as an enzyme producing ligands activating AHR.
#################################################

library(purrr)
library(ggpubr)
library(survival)
library(survminer)
library(RColorBrewer)
library(maxstat)

## Source the functions and parameters files
source("../functions_and_parameters.R")

## Loading the log2 CPM counts data, prior =1 of all TCGA tumors
TCGA_counts <- readRDS("../Zenodo_download/TCGA_counts.rds")
TCGA_counts <- TCGA_counts[tcga_names]
names(TCGA_counts) <- tcga_names

## Loading the GSVA scores of the AHR signature
TCGA_GSVA <- readRDS("../Zenodo_download/TCGA_GSVA_scores_safely.rds")
TCGA_GSVA <- TCGA_GSVA[tcga_names] %>% map(.,function(x){t(x$result)})

## Trp_enzymes
trp_enz <- c("TDO2","IDO1","IL4I1")%>%.[order(.)]

## IL4I1+IDO1+TDO2 Expression levels
Imm_IL4I1_cnts <- map(TCGA_counts, function(x){
  df <- x[rownames(x)%in%trp_enz,]%>%
    t(.)%>%.[,match(trp_enz, colnames(.))]
  data.frame(df) %>% .[order(.$IL4I1),]
})

Imm_IL4I1_cnts <- map2(Imm_IL4I1_cnts, TCGA_GSVA, function(a1,a2){
  idcs <- match(rownames(a1), rownames(a2))
  df <- data.frame(a1, AHR_score=a2[idcs,], stringsAsFactors = F)
})

## Correlating IL4I1 with survival data
# Load TCGA curated clinical data (The 2019 cell reference)
clinical_data <- openxlsx::read.xlsx("../Resources/TCGA_ALL_clinical.xlsx", sheet = 1)

# Create dataframes with Pt_ID identifiers
patient_dfs <- map(Imm_IL4I1_cnts, function(a1){
  df <- data.frame(Pt_ID = rownames(a1)%>%substring(.,1,12) %>% as.character(),a1)
  df
})

# Add the clinical data to each tumor type using the Pt_ID identifier for matching
patient_dfs_2 <- map(patient_dfs, function(a1){
  idx <- match(a1$Pt_ID, clinical_data$bcr_patient_barcode)
  cdata <- clinical_data[idx,c("bcr_patient_barcode","age_at_initial_pathologic_diagnosis","gender",
                               "treatment_outcome_first_course",
                               "OS","OS.time", "DSS", "DSS.time", "DFI","DFI.time","PFI", "PFI.time")]
  colnames(cdata) <- c("bcr_patient_barcode","age","sex","ttt_outcome",
                       "OS","OS.time", "DSS", "DSS.time", "DFI","DFI.time","PFI", "PFI.time")
  cdata
})

# Merge the clinical dataframes and the expression levels
patients_info <- map2(patient_dfs, patient_dfs_2,function(a1,a2){
  df <- data.frame(a1,a2,stringsAsFactors = F)
})

# Add the tumor names to each dataframe
patients_info <- map2(patients_info, names(patients_info),function(a1,a2){
  df <- data.frame(tumor=a2,a1,stringsAsFactors = F)
})

#### Overall survival ####
## Create survival objects based on OS and group IL4I1 using surv_cutpoint
Patients_surv_dfs <- map(patients_info, function(df){
  df$surv.obj <- with(df, Surv(OS.time, OS == 1))
  ct <- surv_cutpoint(df, time="OS.time", event="OS", variables = "IL4I1")
  df$group <- surv_categorize(ct, variables = NULL, labels = c("low", "high"))$IL4I1
  df$group <- factor(df$group,levels = c("low","high"))
  df
})
saveRDS(Patients_surv_dfs, "../Results/RDS/Patients_surv_dfs_IL4I1.rds")

## Cox models of patients split based on IL4I1 expression surv_cutpoint age adjusted
patients_cox_models_IL4I1 <- map(Patients_surv_dfs, function(dat){
  a <- coxph(surv.obj~group+age, data = dat)
})
saveRDS(patients_cox_models_IL4I1, "../Results/RDS/Patients_cox_models_IL4I1_age.rds")

test_coxph_assumption <- map(patients_cox_models_IL4I1, cox.zph)
saveRDS(test_coxph_assumption, "../Results/RDS/test_coxph_assumption_IL4I1.rds")

KM_plots_IL4I1_age <- map(Patients_surv_dfs,function(a){
  df=a
  res.cox <- coxph(surv.obj ~ group + age, data = df)
  p_val_max <- maxstat::maxstat.test(survival::Surv(OS.time, OS) ~ IL4I1, data = df, smethod = "LogRank", pmethod = "Lau92")
  # Create the new data  
  new_df <- with(df, data.frame(group = levels(df$group),
                                age = unlist(map(levels(df$group),function(gl){mean(df$age[df$group==gl],na.rm = T)}))))
  print(new_df)
  fit <- survfit(res.cox, newdata = new_df, data=df)
  ncurve <- ncol(fit$surv)
  ssurv <- tidyr::gather(as.data.frame(fit$surv), "strata", "surv")
  ssurv$time <- rep(fit$time, ncurve )
  p <- ggplot(data = ssurv) + geom_step(aes(time, surv, color = strata))+ggtitle(df$tumor[1])+
    theme_bw()+theme(panel.grid = element_blank())+labs(color="group")+
    scale_color_discrete(name="IL4I1 groups", labels=levels(df$group))+
    annotate(x=200,y=0.5,geom="text",label=paste("p = ",round(p_val_max$p.value,6), sep = ""))
})

pdf("../Results/Figures/IL4I1_OS_surv_km_age_adjusted.pdf", height = 8, width = 12)
KM_plots_IL4I1_age
dev.off()

# Save the IL4I1 cox model representations
patients_cox_models_IL4I1_summaries <- map2(patients_cox_models_IL4I1, tcga_names, function(x,y){
  df <- summary(x)$coefficients
  df <- data.frame(tumor=gsub("TCGA_","",y),vars=rownames(df), df,stringsAsFactors = F)
}) %>% do.call(rbind,.)
patients_cox_models_IL4I1_summaries$vars[patients_cox_models_IL4I1_summaries$vars=="grouphigh"] <- "IL4I1"
patients_cox_models_IL4I1_summaries <- patients_cox_models_IL4I1_summaries[patients_cox_models_IL4I1_summaries$vars=="IL4I1",]
write.table(patients_cox_models_IL4I1_summaries, "../Results/Tables/Patients_cox_models_IL4I1_summaries.txt", sep = "\t", row.names = F)

## Cox models of patients split based on TDO2 expression surv_cutpoint age adjusted
Patients_surv_dfs_TDO2 <- map(patients_info, function(df){
  df$surv.obj <- with(df, Surv(OS.time, OS == 1))
  ct <- surv_cutpoint(df, time="OS.time", event="OS", variables = "TDO2")
  df$group <- surv_categorize(ct, variables = NULL, labels = c("low", "high"))$TDO2
  df$group <- factor(df$group,levels = c("low","high"))
  df
})
saveRDS(Patients_surv_dfs_TDO2, "../Results/RDS/Patients_surv_dfs_TDO2.rds")

patients_cox_models_TDO2 <- map(Patients_surv_dfs_TDO2, function(dat){
  a <- coxph(surv.obj~group+age, data = dat)
})
saveRDS(patients_cox_models_TDO2, "../Results/RDS/Patients_cox_models_TDO2_age.rds")

test_coxph_assumption_TDO2 <- map(patients_cox_models_TDO2, cox.zph)
saveRDS(test_coxph_assumption_TDO2, "../Results/RDS/test_coxph_assumption_TDO2.rds")

KM_plots_TDO2_age <- map(Patients_surv_dfs_TDO2,function(a){
  df=a
  res.cox <- coxph(surv.obj ~ group + age, data = df)
  p_val_max <- maxstat::maxstat.test(survival::Surv(OS.time, OS) ~ TDO2, data = df, smethod = "LogRank", pmethod = "Lau92")
  # Create the new data  
  new_df <- with(df, data.frame(group = levels(df$group),
                                age = unlist(map(levels(df$group),function(gl){mean(df$age[df$group==gl],na.rm = T)}))))
  print(new_df)
  fit <- survfit(res.cox, newdata = new_df, data=df)
  ncurve <- ncol(fit$surv)
  ssurv <- tidyr::gather(as.data.frame(fit$surv), "strata", "surv")
  ssurv$time <- rep(fit$time, ncurve )
  p <- ggplot(data = ssurv) + geom_step(aes(time, surv, color = strata))+ggtitle(df$tumor[1])+
    theme_bw()+theme(panel.grid = element_blank())+labs(color="group")+
    scale_color_discrete(name="IL4I1 groups", labels=levels(df$group))+
    annotate(x=200,y=0.5,geom="text",label=paste("p = ",round(p_val_max$p.value,6),sep = ""))
})

pdf("../Results/Figures/TDO2_OS_surv_km_age_adjusted.pdf", height = 8, width = 12)
KM_plots_TDO2_age
dev.off()

# Save the TDO2 cox model representations
patients_cox_models_TDO2_summaries <- map2(patients_cox_models_TDO2, tcga_names, function(x,y){
  df <- summary(x)$coefficients
  df <- data.frame(tumor=gsub("TCGA_","",y),vars=rownames(df), df, stringsAsFactors = F)
}) %>% do.call(rbind,.)
patients_cox_models_TDO2_summaries$vars[patients_cox_models_TDO2_summaries$vars=="grouphigh"] <- "TDO2"
patients_cox_models_TDO2_summaries <- patients_cox_models_TDO2_summaries[patients_cox_models_TDO2_summaries$vars=="TDO2",]
write.table(patients_cox_models_TDO2_summaries, "../Results/Tables/Patients_cox_models_TDO2_summaries.txt", sep = "\t", row.names = F)

## Cox models of patients split based on IDO1 expression surv_cutpoint age adjusted
Patients_surv_dfs_IDO1 <- map(patients_info, function(df){
  df$surv.obj <- with(df, Surv(OS.time, OS == 1))
  ct <- surv_cutpoint(df, time="OS.time", event="OS", variables = "IDO1")
  df$group <- surv_categorize(ct, variables = NULL, labels = c("low", "high"))$IDO1
  df$group <- factor(df$group,levels = c("low","high"))
  df
})
saveRDS(Patients_surv_dfs_IDO1, "../Results/RDS/Patients_surv_dfs_IDO1.rds")

patients_cox_models_IDO1 <- map(Patients_surv_dfs_IDO1, function(dat){
  a <- coxph(surv.obj~group+age, data = dat)
})
saveRDS(patients_cox_models_IDO1, "../Results/RDS/Patients_cox_models_IDO1_age.rds")

test_coxph_assumption_IDO1 <- map(patients_cox_models_IDO1, cox.zph)
saveRDS(test_coxph_assumption_IDO1, "../Results/RDS/test_coxph_assumption_IDO1.rds")

KM_plots_IDO1_age <- map(Patients_surv_dfs_IDO1,function(a){
  df=a
  res.cox <- coxph(surv.obj ~ group + age, data = df)
  p_val_max <- maxstat::maxstat.test(survival::Surv(OS.time, OS) ~ IDO1, data = df, smethod = "LogRank", pmethod = "Lau92")
  # Create the new data  
  new_df <- with(df, data.frame(group = levels(df$group),
                                age = unlist(map(levels(df$group),function(gl){mean(df$age[df$group==gl],na.rm = T)}))))
  print(new_df)
  fit <- survfit(res.cox, newdata = new_df, data=df)
  ncurve <- ncol(fit$surv)
  ssurv <- tidyr::gather(as.data.frame(fit$surv), "strata", "surv")
  ssurv$time <- rep(fit$time, ncurve )
  p <- ggplot(data = ssurv) + geom_step(aes(time, surv, color = strata))+ggtitle(df$tumor[1])+
    theme_bw()+theme(panel.grid = element_blank())+labs(color="group")+
    scale_color_discrete(name="IL4I1 groups", labels=levels(df$group))+
    annotate(x=200,y=0.5,geom="text",label=paste("p = ",round(p_val_max$p.value,6),sep = ""))
})

pdf("../Results/Figures/IDO1_OS_surv_km_age_adjusted.pdf", height = 8, width = 12)
KM_plots_IDO1_age
dev.off()

# Save the IDO1 cox model representations
patients_cox_models_IDO1_summaries <- map2(patients_cox_models_IDO1, tcga_names, function(x,y){
  df <- summary(x)$coefficients
  df <- data.frame(tumor=gsub("TCGA_","",y),vars=rownames(df), df, stringsAsFactors = F)
}) %>% do.call(rbind,.)
patients_cox_models_IDO1_summaries$vars[patients_cox_models_IDO1_summaries$vars=="grouphigh"] <- "IDO1"
patients_cox_models_IDO1_summaries <- patients_cox_models_IDO1_summaries[patients_cox_models_IDO1_summaries$vars=="IDO1",]
write.table(patients_cox_models_IDO1_summaries, "../Results/Tables/Patients_cox_models_IDO1_summaries.txt", sep = "\t", row.names = F)

## Adding all three dataframes into one summary of all coxph models
coxph_summary_IL4I1_TDO2_IDO1 <- dplyr::bind_rows(patients_cox_models_IL4I1_summaries,patients_cox_models_TDO2_summaries,patients_cox_models_IDO1_summaries)
write.table(coxph_summary_IL4I1_TDO2_IDO1, "../Results/Tables/coxph_summary_IL4I1_TDO2_IDO1.txt", sep = "\t", row.names = F)

## Boxplots of IL4I1 groups and corresponding AHR scores
IL4I1_AHR_OS_groups <- map2(Patients_surv_dfs, tcga_names, function(a1,a2){
  df <- reshape::melt(a1)
  df <- df[df$variable %in% c("IL4I1","AHR_score"),]
  ggboxplot(df, x="group", y="value", fill="group", facet.by = c("variable"),palette = "Set2",title = a2,scales="free", ncol=4)+stat_compare_means()
})

pdf("../Results/Figures/IL4I1_Inhib_mols_OS2.pdf", height = 8, width = 12)
IL4I1_AHR_OS_groups
dev.off()
