#!/bin/bash



base=${1}

source ${base}


sed -i '/^$/d' ${r_log_jobs}PHASE_2_jobs.txt

while (( $(squeue -u adanguy | grep -f ${r_log_jobs}PHASE_2_jobs.txt | wc -l) >= 1 )) 
do    
    sleep 1s
done

# Inclus aussi les mqs de csre
cat ${r_intervals_mqs}*_mqs.txt | sort | uniq > ${r_intervals}mqs.txt


Rscript ${r_scripts}intervals.R ${r_intervals}mqs.txt ${r_amont}SNP_positions.txt ${r_intervals}intervals.txt



 
