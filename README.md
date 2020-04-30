# AHR IL4I1 project scripts
> A collection of scripts used to generate a transcriptional signature for detecting the activation or inhibition states of the Aryl Hydrocarbon Receptor (AHR). In addition, the scripts allow identifying the tryptophan degrading enzyme Interleukin-4 Induced 1 (IL4I1) to activate the AHR through its tryptophan degradation products. These findings and corresponding figures are described in Sadik et al. 2020, "The L-amino acid oxidase IL4I1 is the key enzyme of a novel tryptophan catabolizing pathway and elicits tumor-promoting outcomes through the AHR".

[Introduction]

The Abstract of the paper.

Within this repository we present the (R-) scripts to process and analyse the data as described in Sadik et al. 2020.

## Installation/Prerequisites

All analyses described here were carried out in R and hence can be performed on any operating system (windows, mac or linux). However, running the analysis on a linux distribution is preferred (e.g. ubuntu or centos) because some of the analysis require using multiple cores and are RAM intensive, which would require a computing cluster.

The analyses can be performed in R starting from version 3.3.1 (Bug in Your Hair) and the corresponding bioconductor release. 

libraries

```
R
> install.package("tidyverse")
> install.package("....")

```

## Required files/resources

Describe any additional files that are required that need to be downloaded

@ahmedsadik - perhaps larger files can be linked to zenodo.org

## Paramters to consider

List the scripts and paramters which should be updated to make the workflow run.

@ahmedsadik - perhaps all user definable parameters should be added to a "userParameters.R" file, which is sourced as required?

## Workflow

A few motivating and useful examples of how your product can be used. Spice this up with code blocks and potentially more screenshots.

_For more examples and usage_

| script number | script name | description | 
| --- | --- | --- |
| 0 | functions_and_paramters.R | A script containing functions and parameters pertinent towards more than one other script |
| 1 | TCGA_data.R | The script 'TCGA_data.R' processes [input data] by …. and produces [output data] | 
| 2 | AHR_signature_development.R | The script 'AHR_signature_development.R' processes [input data] by …. and produces [output data] | 
| 3 | AHR_signature_ontologies.R | The script 'AHR_signature_ontologies.R' processes [input data] by …. and produces [output data] | 
| 4 | AHR_signature_validation.R | The script 'AHR_signature_validation.R' processes [input data] by …. and produces [output data] | 
| 5 | GSVA_score.R | The script 'GSVA_score.R' processes [input data] by …. and produces [output data] | 
| 6 | AHR_signature_properties.R | The script 'AHR_signature_properties.R' processes [input data] by …. and produces [output data] | 
| 7 | IDO1_TDO2_associate_w_AHR_signature_activation.R | The script 'IDO1_TDO2_associate_w_AHR_signature_activation.R' processes [input data] by …. and produces [output data] | 
| 8 | TCGA_WGCNA.R | The script 'TCGA_WGCNA.R' processes [input data] by …. and produces [output data] | 
| 9 | WGCNA_modules_Associating_w_AHR_activation.R | The script 'WGCNA_modules_Associating_w_AHR_activation.R' processes [input data] by …. and produces [output data] | 
| 10 | WGCNA_modules_Associating_w_AHR_activation_Gene_ontologies.R | The script 'WGCNA_modules_Associating_w_AHR_activation_Gene_ontologies.R' processes [input data] by …. and produces [output data] | 
| 11 | IDO1_TDO2_circos_WGCNA_modules_associating_w_AHR_signature_activation.R | The script 'IDO1_TDO2_circos_WGCNA_modules_associating_w_AHR_signature_activation.R' processes [input data] by …. and produces [output data] | 
| 12 | IL4I1_activates_AHR.R | The script 'IL4I1_activates_AHR.R' processes [input data] by …. and produces [output data] | 
| 13 | GBM_WGCNA_modules_and_selected_TRP_enzymes_circos_plot.R | The script 'GBM_WGCNA_modules_and_selected_TRP_enzymes_circos_plot.R' processes [input data] by …. and produces [output data] | 
| 14 | WGCNA_modules_Consensus_Clustering.R | The script 'WGCNA_modules_Consensus_Clustering.R' processes [input data] by …. and produces [output data] | 
| 15 | Validating_cluster_separation_using_PCA.R | The script 'Validating_cluster_separation_using_PCA.R' processes [input data] by …. and produces [output data] | 
| 16 | Differential_Gene_Expression_of CC_groups.R | The script 'Differential_Gene_Expression_of CC_groups.R' processes [input data] by …. and produces [output data] | 
| 17 | ROAST_GSA.R | The script 'ROAST_GSA.R' processes [input data] by …. and produces [output data] | 
| 18 | GSVA_boxplots_selected_cluster_assignments.R | The script 'GSVA_boxplots_selected_cluster_assignments.R' processes [input data] by …. and produces [output data] | 
| 19 | IL4I1_levels_in_AHR_subtypes.R | The script 'IL4I1_levels_in_AHR_subtypes.R' processes [input data] by …. and produces [output data] | 
| 20 | Heatmaps_modules_posModulesGenes.R | The script 'Heatmaps_modules_posModulesGenes.R' processes [input data] by …. and produces [output data] | 
| 21 | Infilteration_gene_signatures.R | The script 'Infilteration_gene_signatures.R' processes [input data] by …. and produces [output data] | 
| 22 | GSVA_scores_Immune_cell_types.R | The script 'GSVA_scores_Immune_cell_types.R' processes [input data] by …. and produces [output data] | 
| 23 | Immune_infiltarting_populations_distribution_AHR_clusters.R | The script 'Immune_infiltarting_populations_distribution_AHR_clusters.R' processes [input data] by …. and produces [output data] | 
| 24 | Immune_infiltration_per_patient_heatmaps.R | The script 'Immune_infiltration_per_patient_heatmaps.R' processes [input data] by …. and produces [output data] | 
| 25 | Initial_SKCM_script.R | The script 'Initial_SKCM_script.R' processes [input data] by …. and produces [output data] | 

## The workflow

0) Start by cloning the reposistory and external data from zenodo.org
```
mkdir -p /path/to/AHR_workflow
cd /path/to/AHR_workflow
git clone https://github.com/ahmedasadik/Project_AHR_LAAO.git
mkdir zenodo_data
cd zenodo_data
wget https://zenodo.org/record/[doi-goes-here]/files/[file1.txt]
wget https://zenodo.org/record/[doi-goes-here]/files/[file2.txt]
cd ..
```

Make subdirectories required for the R script outputs
```
cd Project_AHR_LAAO
mkdir Figures
mkdir IDO1
mkdir in
mkdir RDS
mkdir sckm_metastasis
mkdir signed
mkdir SKCM
mkdir ...
cd ..
```

Make sure that the correct version of R is installed and download expected libraries

```
> R --version
R version 3.3.1 (2016-06-21) -- "Bug in Your Hair"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under the terms of the
GNU General Public License versions 2 or 3.
For more information about these matters see
http://www.gnu.org/licenses/.

> R
install.packages("AffyPipelineHG2ST")
install.packages("biomaRt")
install.packages("circlize")
install.packages("clusterProfiler")
install.packages("colorspace")
install.packages("ComplexHeatmap")
install.packages("ConsensusClusterPlus")
install.packages("corrplot")
install.packages("dplyr")
install.packages("edgeR")
install.packages("factoextra")
install.packages("FactoMineR")
install.packages("GetoptLong")
install.packages("GGally")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("globaltest")
install.packages("gplots")
install.packages("GSVA")
install.packages("miceadds")
install.packages("org.Hs.eg.db")
install.packages("parallel")
install.packages("psych")
install.packages("purrr")
install.packages("RColorBrewer")
install.packages("TCGAbiolinks")
install.packages("WGCNA")
```

1) TCGA_data.R

Here we processes [input data] by …. and produces [output data]

```
> cd Project_AHR_LAAO
> cd AHR_scripts
> R -f TCGA_data.R

...

[code block?]
...
```

At the end we produce XXX, YYY and ZZZ output files. The output of YYY is visualised as:

![](TCGA_data_plot.png)

2) ...


## Authors, copyright, license

Dr Ahmed Sadik –  a.sadik@dkfz.de - [https://github.com/ahmedsadik](https://github.com/ahmedsadik/)

Dr Christiana Opitz - c.opitz@dkfz.de - [https://www.dkfz.de/en/brain-cancer-metabolism/index.php](https://www.dkfz.de/en/brain-cancer-metabolism/index.php)

[Copyright statement]

Distributed under the XYZ license. See ``LICENSE`` for more information.

## Session info

```
R version 3.3.1 (2016-06-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: openSUSE 13.1 (Bottle) (x86_64)

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base
```
