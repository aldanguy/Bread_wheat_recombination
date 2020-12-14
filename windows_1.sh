#!/bin/bash


# Define windows for PHASE inputs

base=${1}


source ${base}

rm ${r_PHASE_windows_description}*
rm ${r_PHASE_windows_list}*


for p in ${pop[*]}
do

    for c in ${chr[*]}
    do
    
    echo ${p}
    
    echo ${c}


    job_out=${r_log_windows}w_${p}_${c}.out

    job_name=${p}_${c}

    job=$(sbatch -o ${job_out} -J ${job_name} --time=00:30:00 --mem=10G --parsable ${r_scripts}windows_2.sh ${base} ${p} ${c})
    
    
    echo "${job}" >> ${r_log_jobs}windows_2_jobs.txt
    
    done
    
done





sed -i '/^$/d' ${r_log_jobs}windows_2_jobs.txt
while (( $(squeue -u adanguy | grep -f ${r_log_jobs}windows_2_jobs.txt | wc -l) >= 1 )) 
do    
    sleep 1s
done






