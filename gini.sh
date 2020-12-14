#!/bin/bash

base=${1}


source ${base}



Rscript ${r_scripts}gini_populations.R ${r_PHASE}PHASE_summary_outputs.txt ${r_intervals}intervals.txt ${r_sources}csre_genetic_map.txt ${r_gini}pour_gini.txt ${r_graphs}gini_populations.png ${r_gini}gini_populations.txt

Rscript ${r_scripts}gini_genetic_distance.R ${r_PHASE}PHASE_summary_outputs.txt ${r_sources}csre_genetic_map.txt ${r_HR}HR.txt ${r_graphs}gini_genome.png ${r_gini}gini_distance_genetique.txt

