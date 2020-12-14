#!/bin/bash

## For each SNP of the dataset, provide an ID, a physical position and a origin

base=${1}
source ${base}

Rscript ${r_scripts}SNP_positions.R ${pos} ${SNP_category} ${chr_names} ${chr_tab} ${csre_others_SNP} ${PLINK_initial}.map ${r_sources}SNP_positions.txt

