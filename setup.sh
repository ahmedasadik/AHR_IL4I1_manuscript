#!/bin/bash

mkdir -p $HOME/AHR_IL4I1_workflows
cd $HOME/AHR_IL4I1_workflows
git clone https://github.com/ahmedasadik/AHR_IL4I1_manuscript.git

cd AHR_IL4I1_manuscript

# creating download directories
mkdir Zenodo_download
# please download the files using the shared link that you received till it becomes public
# cd Zenodo_download
# download necessary files from the zendo repository
# wget https://zenodo.org/record/[doi-goes-here]/files/[file1.txt]
# wget https://zenodo.org/record/[doi-goes-here]/files/[file2.txt]
# wget https://zenodo.org/record/[doi-goes-here]/files/[file3.txt]
# cd ..

# creating workflow directories
mkdir Results TCGAbiolinks_downloads
mkdir ./Results/AHR_signature ./Results/AHR_signature_validations ./Results/RDS ./Results/GlobalTest ./Results/TCGA_Circos_plots ./Results/TCGA_Circos_plots/Circos_plots_2 ./Results/TCGA_Circos_plots/Circos_plots_7 ./Results/Tables ./Results/Figures ./Results/GO_AAMs ./Results/GO_AAMs/BP_sim ./Results/GO_AAMs/BP
mkdir ./TCGAbiolinks_downloads/TCGA_counts ./TCGAbiolinks_downloads/TCGA_counts/RDats ./TCGAbiolinks_downloads/TCGA_FPKM ./TCGAbiolinks_downloads/TCGA_FPKM/RDats
