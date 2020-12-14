#!/bin/bash

# Estimate genetic positions of SNP wich are unmpapped in CsRe population, based on their physical positions and the physical and genetic positions of mapped CsRe SNP

base=${1}
source ${base}

Rscript ${r_scripts}interpolation.R ${r_sources}csre_genetic_map.txt ${r_sources}SNP_positions.txt ${r_graphs_interpolation} ${r_PLINK}update_posSNP.txt ${r_PLINK}update_genSNP.txt ${r_PLINK}SNP_PLINK.txt ${r_graphs}bayesian_and_interpolation.png
