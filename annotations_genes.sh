#!/bin/bash

base=${1}


source ${base}



cat ${r_amont}iwgsc_refseqv1.0_HighConf_UTR_2017May05.gff3 | grep -v "##" > ${r_amont}annotation.txt


cut -f9 ${r_amont}annotation.txt | cut -f1 -d";" | sed "s|ID=||" | cut -f1 -d"." > ${r_amont}tempgenes.txt


# cat ${r_amont}annotation.txt | cut -f1-8 | paste -d ${r_amont}tempgenes.txt

paste ${r_amont}annotation.txt ${r_amont}tempgenes.txt | cut -f1,3,4,5,10 > ${r_amont}tempgenes2.txt

sed "s|chr||" ${r_amont}tempgenes2.txt  > ${r_amont}tempgenes.txt
cp ${r_amont}tempgenes.txt ${r_amont}annotation.txt
rm ${r_amont}tempgenes2.txt
rm ${r_amont}tempgenes.txt


Rscript ${r_scripts}annotations_genes.R ${r_amont}annotation.txt ${r_HR}HR.txt ${chr_tab} ${r_PHASE}PHASE_summary_outputs.txt ${r_HR}annotations_genes_HR.txt ${r_graphs}annotation_genes_HR.png
