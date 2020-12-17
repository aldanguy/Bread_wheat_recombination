#!/bin/bash

## Bayesian genetic map for CsRe

base=${1}
source ${base}



Rscript ${r_scripts}csre_bayesian.R ${r_csre}csre_genetic_map.txt ${r_csre}csre_crossovers_sampling.txt ${r_graphs} ${r_graphs}csre_recombination_rates_ECDF.png ${r_csre}csre_mqs.txt


