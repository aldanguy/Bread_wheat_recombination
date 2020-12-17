#!/bin/bash

## Passeport data about landraces : ID, STRUCTURE results

base=${1}
source ${base}


Rscript ${r_scripts}landraces.R ${accessions} ${structure} ${r_landraces}landraces.txt ${K4} ${K8}

# extract "LINE" id

cat ${r_landraces}landraces.txt | cut -f4 | sed "s|^|1\t|" | tail -n +2 > ${r_PLINK}landraces_PLINK.txt
