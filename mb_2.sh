#!/bin/bash

base=${1}


source ${base}

p=${2}

c=${3}

carte=${4}

titre=${5}




Rscript ${r_scripts}mb_landraces.R ${r_PHASE}${p}_${c}_PHASE_raw_outputs.txt ${chr_tab} ${carte} ${r_mb}${titre}.txt


