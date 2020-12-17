#!/bin/bash



base=${1}

p=${2}

c=${3}

source ${base}


Rscript ${r_scripts}windows_2.R ${p} ${c} ${r_PLINK} ${r_amont}SNP_positions.txt ${r_PHASE_windows_list} ${r_PHASE_windows_description}w_${p}_${c}.txt



