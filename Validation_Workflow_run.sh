#!/bin/bash

# download the necessary files from the zendo repository and untar it here
# it comes with the right directory tree needed for the validation scripts
# wget https://zenodo.org/record/[doi-goes-here]/files/[file1.txt]
# tar -xf Validation_datasets.tar

## Run the IL4I1_scripts
cd Validation_scripts
./GSE102045_FICZ_Th17.R
./GSE109576_CH22_TCDD_A549.R
./GSE25272_Kyn_U87.R
./GSE28878_TCDD_HepG2.R
./GSE32026_TCDD_hMADs.R
./GSE48843_SR1_AML.R
./GSE67093_SR1_HSC.R
./GSE98515_TCDD_MCF7.R
cd ..