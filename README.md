# AHR IL4I1 project scripts
A collection of scripts used to generate a transcriptional signature for detecting the activation or inhibition states of the Aryl Hydrocarbon Receptor (AHR). In addition, the scripts allow identifying the tryptophan degrading enzyme Interleukin-4 Induced 1 (IL4I1) to activate the AHR through its tryptophan degradation products. These findings and corresponding figures are described in Sadik et al. 2020, "The L-amino acid oxidase IL4I1 is the key enzyme of a novel tryptophan catabolizing pathway and elicits tumor-promoting outcomes through the AHR".

Within this repository we present the (R-) scripts to process and analyse the data as described in Sadik et al. 2020.

## Installation/Prerequisites

All analyses described here were carried out in R and hence can be performed on any operating system (windows, mac or linux). However, running the analysis on a linux distribution is preferred (e.g. ubuntu or centos) because some of the analysis require using multiple cores and are RAM intensive, which would require a computing cluster.

All of the scripts were independently tested on a clean conda software environment to ensure reproducibility.

## Required files/resources

Before you start, please be sure that you have R version 3.6.1 and bioconductor installed on your machine. Then proceed with downloading the setup.sh file and then run it from the command line. This will create a folder in your home directory and will clone the repository, download all required files from zenoodo and create all the necessary output directories. For the microarrays that were used for establishing the relations between the AHR and IL4I1 download the raw cel files and place them in the respective directory (e.g. './Results/IL4I1_microarrays/I3P_SR1/CEL'). 

Before running any workflows, be sure that your current active directory is the ./AHR_IL4I1_manuscript directory that was cloned from github.
To run the different workflows, start by running the Workflow_run.sh file which will recreat all the figures of the analyses establishing the relation between AHR and IL4I1. The validation of the AHR signature results can be generated by running the Validation_Workflow_run.sh file. If you would like to regenrate the TCGA data objects that were downloaded from zenodo, run the generate_TCGA_data.sh, but be aware that this requires using a computing cluster because many of the scripts require very large RAM and CPU power.
