#!/bin/bash

## CONDA SETUP for R.3.6.0/3.6.1 (although @ahmedsadik did local dev on 3.6.3)
# conda create --name AHR_IL4I1_env
# conda activate AHR_IL4I1_env
# conda install r=3.6.0 python mofapy

##  Installing MOFA dependencies
# you can install this via several methods: https://github.com/bioFAM/MOFA#installation
# conda: conda install -c bioconda mofapy
# pip: pip install mofapy

mkdir -p $HOME/AHR_IL4I1_workflows
cd $HOME/AHR_IL4I1_workflows
git clone https://github.com/ahmedasadik/AHR_IL4I1_manuscript.git

cd AHR_IL4I1_manuscript

# creating download directories
mkdir Zenodo_download
# please download the files using the shared link that you received till it becomes public
cd Zenodo_download
# download necessary files from the zendo repository
wget https://zenodo.org/record/3778914/files/GT_GSVA_overlap_MEs.rds
wget https://zenodo.org/record/3778914/files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct
wget https://zenodo.org/record/3778914/files/GTEX_annot
wget https://zenodo.org/record/3778914/files/IPS_cell_types_gene_lists.rds
wget https://zenodo.org/record/3778914/files/TCGA_coldata.rds
wget https://zenodo.org/record/3778914/files/TCGA_counts.rds
wget https://zenodo.org/record/3778914/files/TCGA_DGE.rds
wget https://zenodo.org/record/3778914/files/TCGA_DGE_voom.rds
wget https://zenodo.org/record/3778914/files/TCGA_DGE_voom_annot.rds
wget https://zenodo.org/record/3778914/files/TCGA_fpkms_all.rds
wget https://zenodo.org/record/3778914/files/TCGA_GSVA_scores_safely.rds
wget https://zenodo.org/record/3778914/files/TCGA_GSVA_voom_28_cell_types_IPS.rds
wget https://zenodo.org/record/3778914/files/TCGA_MEs_GSVA_cors.rds
wget https://zenodo.org/record/3778914/files/TCGA_MEs_GSVA_GT_overlap_nodir.rds
wget https://zenodo.org/record/3778914/files/TCGA_MEs_GSVA_GTs_not_null_sig05_nodir.rds
wget https://zenodo.org/record/3778914/files/TCGA_modules_genes.rds
wget https://zenodo.org/record/3778914/files/TCGA_sft.rds
wget https://zenodo.org/record/3778914/files/TCGA_TOMS_bicor.rds
wget https://zenodo.org/record/3778914/files/TCGA_TPMs.rds
cd ..

# creating workflow directories
mkdir Results TCGAbiolinks_downloads
mkdir ./Results/AHR_signature ./Results/RDS ./Results/GlobalTest ./Results/TCGA_Circos_plots ./Results/TCGA_Circos_plots/Circos_plots_2 ./Results/TCGA_Circos_plots/Circos_plots_7 ./Results/Tables ./Results/Figures ./Results/GO_AAMs ./Results/GO_AAMs/BP_sim ./Results/GO_AAMs/BP
mkdir ./Results/Validations ./Results/Validations/GSE102045_Th17 ./Results/Validations/GSE109576_A549 ./Results/Validations/GSE25272_U87 ./Results/Validations/GSE28878_HepG2 ./Results/Validations/GSE32026_hMADs ./Results/Validations/GSE48843_AML ./Results/Validations/GSE67093_HSC ./Results/Validations/GSE98515_MCF7
mkdir ./TCGAbiolinks_downloads/TCGA_counts ./TCGAbiolinks_downloads/TCGA_counts/RDats ./TCGAbiolinks_downloads/TCGA_FPKM ./TCGAbiolinks_downloads/TCGA_FPKM/RDats
mkdir ./Results/IL4I1_microarrays ./Results/IL4I1_microarrays/CAS1 ./Results/IL4I1_microarrays/CAS1/CEL ./Results/IL4I1_microarrays/Metabolites ./Results/IL4I1_microarrays/Metabolites/CEL ./Results/IL4I1_microarrays/I3A ./Results/IL4I1_microarrays/I3A/CEL
mkdir ./Results/IL4I1_microarrays/I3P_SR1 ./Results/IL4I1_microarrays/I3P_SR1/CEL ./Results/IL4I1_microarrays/IL4I1_shAHR ./Results/IL4I1_microarrays/IL4I1_shAHR/CEL ./Results/IL4I1_microarrays/KynA_FICZ ./Results/IL4I1_microarrays/KynA_FICZ/CEL
mkdir ./Results/GSE66858

# install R libraries and dependencies
./install_required_packages.R
