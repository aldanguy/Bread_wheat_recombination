#!/bin/bash


base=${1}

source ${base}


Rscript ${r_scripts}select_window_for_DL_graph_1.R ${r_HR}HR.txt ${r_PHASE}windows.txt ${r_PHASE}window_for_DL.txt

p=WE

cat ${r_PHASE}window_for_DL.txt | while read line 
do

plink --file ${r_PLINK}${p} --extract ${r_PHASE_windows_list}w_${line}.txt --r2 --ld-snp-list ${r_PHASE_windows_list}w_${line}.txt --ld-window-kb 999 --ld-window 99999 --ld-window-r2 0 --noweb --out ${r_PHASE}w_${line}

rm ${r_PHASE}w_${line}.log


done


Rscript ${r_scripts}select_window_for_DL_graph_2.R ${r_PHASE}w_${line}.ld ${r_HR}HR.txt ${r_PHASE}windows.txt ${r_PHASE}PHASE_summary_outputs.txt ${r_graph}LD.tiff
