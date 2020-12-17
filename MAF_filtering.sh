#!/bin/bash

## Extract polymorphic SNP per population : MAF >= 3%

base=${1}

source ${base}

for p in ${pop[*]}
do

    echo ${p}
    


    plink --file ${r_landraces}landraces --recode --noweb --extract ${r_PLINK}SNP_PLINK_2.txt --keep ${r_PLINK}${p}.txt --maf ${rate_maf_min} --out ${r_PLINK}${p} >>${r_log}MAF_filtering.out
  


done

rm ${r_PLINK}*.log


Rscript ${r_scripts}MAF_filtering.R ${r_landraces}landraces.txt ${r_PLINK}WE.map ${r_PLINK}EE.map ${r_PLINK}WA.map ${r_PLINK}EA.map ${r_amont}SNP_positions.txt ${r_PLINK} ${r_PLINK}update_FID.txt ${r_PLINK}SNP_PLINK_3.txt ${r_PLINK}landraces_PLINK_3.txt

plink --file ${r_landraces}landraces --recode --noweb --keep ${r_PLINK}landraces_PLINK_3.txt --extract ${r_PLINK}SNP_PLINK_3.txt --out ${r_landraces}landraces_final >>${r_log}MAF_filtering.out

plink --file ${r_landraces}landraces_final --recode --noweb --update-ids ${r_PLINK}update_FID.txt --out ${r_landraces}landraces_final >>${r_log}MAF_filtering.out

plink --file ${r_landraces}landraces_final --noweb --make-bed --out ${r_landraces}landraces_final >>${r_log}MAF_filtering.out

    
rm ${r_landraces}*.log
