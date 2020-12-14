#!/bin/bash

## CsRe genotyping, crossovers identification, sampling of crossovers positions, empirical genetic map 


base=${1}

source ${base}


Rscript ${r_scripts}csre_genotyping.R ${raw_csre_genotyping} ${csre_others_SNP} ${r_sources}SNP_positions.txt ${r_sources}csre_genotyping.txt


Rscript ${r_scripts}csre_crossovers.R ${r_sources}csre_genotyping.txt ${r_sources}scaffolds.txt ${r_sources}csre_genetic_map.txt ${r_csre}csre_crossovers_sampling.txt ${nb_posterior} ${r_csre}csre_crossovers_positions.txt





