#!/bin/bash


base=${1}

source ${base}



Rscript ${r_scripts}landraces_structuration.R ${r_landraces}landraces.txt ${dis} ${r_PLINK}landraces_PLINK_2.txt ${r_PLINK} ${r_graphs}





