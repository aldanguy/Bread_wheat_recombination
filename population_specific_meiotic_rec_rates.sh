#!/bin/bash

base=${1}


source ${base}


Rscript ${r_scripts}population_specific_meiotic_rec_rates.R ${r_PHASE}PHASE_summary_outputs.txt ${r_csre}csre_genetic_map.txt ${r_maps_pop}landraces_genetic_maps.txt ${r_graphs}proportionnality_constant.png

