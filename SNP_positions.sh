#!/bin/bash

## For each SNP of the dataset, provide an ID, a physical position and a origin

base=${1}
source ${base}


# You can't run this script because some files contain private SNP. The script consists in linking SNP ID and their physical positions. The output needed for futher scripts is SNP_positions.txt
# Rscript ${r_scripts}SNP_positions.R ${r_amont}Vraies_positions_marqueurs.txt ${r_amont}Axiom_TaBW420_OrigineSNPs.csv ${chr_names} ${chr_tab} ${r_amont}CartPolyHigh-Oriented_def_def_def.txt ${PLINK_initial}.map ${r_amont}SNP_positions.txt

