#!/bin/bash

## CsRe genotyping, crossovers identification, sampling of crossovers positions, empirical genetic map 


base=${1}

source ${base}

# You can't run this script because MatriceGenoHead.CorrectedName.UniqSorted.tab and CartPolyHigh-Oriented_def_def_def.txt are not published. The script consist only in filtering some private SNP. The output needed for futher scripts is csre_genotyping.txt
# Rscript ${r_scripts}csre_genotyping.R ${r_amont}MatriceGenoHead.CorrectedName.UniqSorted.tab ${r_amont}CartPolyHigh-Oriented_def_def_def.txt ${r_amont}SNP_positions.txt ${r_amont}csre_genotyping.txt


Rscript ${r_scripts}csre_crossovers.R ${csre_genotyping} ${r_scaffolds}scaffolds.txt ${r_csre}csre_genetic_map.txt ${r_csre}csre_crossovers_sampling.txt ${nb_posterior} ${r_csre}csre_crossovers_positions.txt





