#!/bin/bash


base=${1}

source ${base}



p=${2}

c=${3}

f=${4}



Rscript ${r_scripts}intervals_PHASE_2.R ${r_intervals}intervals.txt ${r_PHASE_output_format}${p}_${c}_${f}_PHASE_raw_outputs.txt ${r_intervals_PHASE}${p}_${c}_${f}_intervals.txt



