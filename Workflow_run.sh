#!/bin/bash

## Run the AHR signature scripts first
cd AHR_signature_scripts
./IDO1_TDO2_associate_w_AHR_signature_activation.R
./WGCNA_modules_Associating_w_AHR_activation_AAMs.R
cd ..

## Run the IL4I1_scripts
cd IL4I1_scripts
./IL4I1_activates_AHR.R
./IL4I1_and_AHR_activity_in_metastatic_and_primary_melanoma.R
./IL4I1_immune_associations.R
./IL4I1_in_TCGA_and_GTEX_tissues.R
./IL4I1_survival_outcome_TCGA.R
./Riaz_etal_ICB_dataset.R
./Weisser_GSE58211_relapsed_CLL_dataset.R
cd ..

## Run the IL4I1_scripts_microarrays
cd IL4I1_scripts_microarrays
./CAS1.R
./HPP_PP_I3P_metabolites.R
./I3C.R
./I3P_SR1.R
./IL4I1_shAHR.R
./KynA_FICZ.R
cd ..