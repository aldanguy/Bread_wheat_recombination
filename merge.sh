#!/bin/bash



base=${1}

source ${base}


sed -i '/^$/d' ${r_log_jobs}intervals_PHASE_2_jobs.txt

while (( $(squeue -u adanguy | grep -f ${r_log_jobs}intervals_PHASE_2_jobs.txt | wc -l) >= 1 )) 
do    
    sleep 1s
done











## sampling in intervals
echo "intervals"

k=0
for f in ${r_intervals_PHASE}*intervals.txt
do 


    if ((k==0 ))
        
    then 
        
        cat ${f} > ${r_intervals}intervals.txt
            # rm ${f}
        
    else
        
        tail -n+2 ${f} >> ${r_intervals}intervals.txt
            # rm ${f}

    fi
        
    k=$((k +1))
        
          
done

Rscript ${r_scripts}merge_PHASE_intervals.R ${r_intervals}intervals.txt ${r_intervals}intervals.txt




## summary outputs of PHASE
echo "summary"
k=0
for f in ${r_PHASE_output_summary}*_PHASE_summary_outputs.txt
do 

    if ((k==0 ))
        
    then 
        
        cat ${f} > ${r_PHASE}PHASE_summary_outputs.txt
            # rm ${f}
        
    else
        
        tail -n+2 ${f} >> ${r_PHASE}PHASE_summary_outputs.txt
            # rm ${f}

    fi
        
    k=$((k +1))
        
          
done



## windows
echo "w"
k=0
for f in ${r_PHASE_windows_description}w*.txt
do 

    if ((k==0 ))
        
    then 
        
        cat ${f} > ${r_PHASE}windows.txt
            # rm ${f}
        
    else
        
        tail -n+2 ${f} >> ${r_PHASE}windows.txt
            # rm ${f}

    fi
        
    k=$((k +1))
        
          
done



## PHASE


echo "time"

cat ${r_PHASE_time}* > ${r_PHASE}time.txt

Rscript ${r_scripts}merge_PHASE_times.R ${r_PHASE}windows.txt ${r_PHASE}time.txt

rm ${r_PHASE}time.txt




## PHASE raw outputs

echo "raw"
k=0
for f in ${r_PHASE_output_format}*_PHASE_raw_outputs.txt 
do 

    if ((k==0 ))
        
    then 
        
        cat ${f} > ${r_PHASE}PHASE_raw_outputs.txt 
            # rm ${f}
        
    else
        
        tail -n+2 ${f} >> ${r_PHASE}PHASE_raw_outputs.txt 
            # rm ${f}

    fi
        
    k=$((k +1))
        
          
done



for c in ${chr[*]}
do

    echo ${c}

    for p in ${pop[*]}
    do
    
    echo ${p}
    
    awk -v c="${c}" -v p="${p}" 'BEGIN {FS="\t"} (NR==1) {print $0} ($2==c && $1==p) {print $0}' ${r_PHASE}PHASE_raw_outputs.txt  > ${r_PHASE}${p}_${c}_PHASE_raw_outputs.txt 
    
done

done

rm ${r_PHASE}PHASE_raw_outputs.txt

