#!/bin/bash

base=${1}


source ${base}


ID=${3}


marges=${2}


echo ${ID}

Rscript ${r_scripts}HR_permutation.R ${r_PHASE}PHASE_summary_outputs.txt ${r_HR}HR.txt ${chr_tab} ${ID} ${r_HR_permutations}permutation_${ID}.txt 


Rscript ${r_scripts}HR_coloc.R permutation ${r_HR_permutations}permutation_${ID}.txt ${r_PHASE}PHASE_summary_outputs.txt ${r_HR_permutations}permutation_${ID}.txt ${marges}


