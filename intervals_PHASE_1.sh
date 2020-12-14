#!/bin/bash



base=${1}

source ${base}

rm ${r_intervals_PHASE}*


for c in ${chr[*]}
do


    for p in ${pop[*]}
    do
    
    f=1

    


   
 max_fen=$(tail -n1 ${r_PHASE_windows_description}w_${p}_${c}.txt | cut -f9 | sed "s/^.*_//g")   
    
while (( ${f} <= ${max_fen} ))
do
	
	nb_jobs_squeue=$(squeue -u adanguy | wc -l)

	
	
	if (( ${nb_jobs_squeue} <= ${nb_jobs_allowed} )) 
	
	then
    
    
  
	echo ${p}_${c}_${f}
    
        
    job_out=${r_log_intervals}${p}_${c}_${f}_intervals.out
    
    
    job_name=${f}_${c}_${p}
	
    
    job=$(sbatch -o ${job_out} -J ${job_name} --mem=20G --parsable ${r_scripts}intervals_PHASE_2.sh ${base} ${p} ${c} ${f})
    
    
    echo "${job}" >> ${r_log_jobs}intervals_PHASE_2_jobs.txt
		
		
    f=$((f +1))
		
	else
	
		sleep 1s
		
	fi
    
    
done

done


done


sed -i '/^$/d' ${r_log_jobs}intervals_PHASE_2_jobs.txt
while (( $(squeue -u adanguy | grep -f ${r_log_jobs}intervals_PHASE_2_jobs.txt | wc -l) >= 1 )) 
do    
    sleep 1s
done
 



rep_out=${r_log_intervals}
rep_sortie=${r_intervals_PHASE}
# Recherche de pb
beug=$(grep -Ril -e "segfault" -e "error" ${rep_out})  
if [ ! -z "${beug}" ]
then
beug2=$(echo ${beug} | sed "s|${rep_out}||g" | sed "s|.out||g")
# supression fichiers sorties (rq : le fichier log et sortie doivent s'appeler pareil sauf les repertoires et extensions)
supression=$(echo ${beug} | sed "s|${rep_out}|${rep_sortie}|g" | sed "s|.out|.txt|g")
echo ${beug2}
rm ${supression}
fi
