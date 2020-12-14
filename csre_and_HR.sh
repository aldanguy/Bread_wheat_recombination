#!/bin/bash

base=${1}


source ${base}


Rscript ${r_scripts}csre_and_HR.R ${r_HR}HR.txt ${r_csre}map_4mb_csre.txt ${r_csre}csre_crossovers_positions.txt ${r_graphs}CsRe_HR.png



map_4mb_csre=${r_csre}map_4mb_csre.txt
csre_crossovers_positions=${r_csre}csre_crossovers_positions.txt


source ${base_common_SNP}



while [ ! -f  ${r_HR}HR.txt ] ; 
do 
    sleep 1 m
    
done


Rscript ${r_scripts}csre_and_HR.R ${r_HR}HR.txt ${map_4mb_csre} ${csre_crossovers_positions} ${r_graphs}CsRe_HR.png
