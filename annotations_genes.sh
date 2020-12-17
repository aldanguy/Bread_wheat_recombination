#!/bin/bash

base=${1}


source ${base}



cat ${amont_annotation} | grep -v "##" > ${r_annotation}annotation.txt


cut -f9 ${r_annotation}annotation.txt | cut -f1 -d";" | sed "s|ID=||" | cut -f1 -d"." > ${r_annotation}tempgenes.txt


# cat ${r_annotation}annotation.txt | cut -f1-8 | paste -d ${r_annotation}tempgenes.txt

paste ${r_annotation}annotation.txt ${r_annotation}tempgenes.txt | cut -f1,3,4,5,10 > ${r_annotation}tempgenes2.txt

sed "s|chr||" ${r_annotation}tempgenes2.txt  > ${r_annotation}tempgenes.txt
cp ${r_annotation}tempgenes.txt ${r_annotation}annotation.txt
rm ${r_annotation}tempgenes2.txt
rm ${r_annotation}tempgenes.txt


Rscript ${r_scripts}annotations_genes.R ${r_annotation}annotation.txt ${r_HR}HR.txt ${chr_tab} ${r_PHASE}PHASE_summary_outputs.txt ${r_HR}annotations_genes_HR.txt ${r_graphs}annotation_genes_HR.png
