#!/bin/bash


## Extract genotyping data : filtered markers and filtered landraces (het : 5% max ; NA : 10% max ; polymorphic markers)


base=${1}



source ${base}

pop=landraces


plink --file ${PLINK_initial} --recode --keep ${r_PLINK}landraces_PLINK.txt --extract ${r_PLINK}SNP_PLINK.txt --update-map ${r_PLINK}update_posSNP.txt --noweb --out ${r_landraces}landraces >>${r_log}PLINK.out

plink --file ${r_landraces}landraces --recode --update-map ${r_PLINK}update_genSNP.txt --update-cm --noweb --out ${r_landraces}landraces >>${r_log}PLINK.out

plink --file ${r_landraces}landraces --recode --maf 0.000000001 --noweb --geno ${rate_NA_max} --mind ${rate_NA_max} --out ${r_landraces}landraces >>${r_log}PLINK.out


plink --file ${r_landraces}landraces --recode --freq --noweb --out ${r_landraces}landraces >>${r_log}PLINK.out




Rscript ${r_scripts}rate_heterozygotie_NA.R ${pop} ${r_landraces}landraces.ped ${r_landraces}landraces.map ${r_landraces}landraces.frq ${r_amont}SNP_positions.txt ${r_PLINK}SNP_PLINK_2.txt ${r_PLINK}landraces_PLINK_2.txt ${r_PLINK}NA.txt ${r_PLINK}NA.txt ${rate_het_max} 1 0



rm ${r_landraces}landraces.log
rm ${r_PLINK}NA.txt


