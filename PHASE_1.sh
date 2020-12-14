#!/bin/bash


base=${1}

source ${base}

<<COMMENTS

rm ${r_log_PHASE}PHASE_jobs.txt
rm ${r_log_PHASE}*
rm ${r_PHASE_input_inp}*
rm ${r_PHASE_time}*
rm ${r_PHASE_output_raw}*
rm ${r_PHASE_output_format}*
rm ${r_PHASE_output_summary}*

COMMENTS


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
  
	job_out=${r_log_PHASE}PHASE_${p}_${c}_${f}.out
    
    job_name=${f}_${c}_${p}
	    
    job_PHASE=$(sbatch -o ${job_out} -J ${job_name} --mem=20G --parsable ${r_scripts}PHASE_2.sh ${base} ${p} ${c} ${f}) # change for unlimitq
    
    
    echo "${job_PHASE}" >> ${r_log_jobs}PHASE_2_jobs.txt

		
    f=$((${f} +1))
		
	else
	
		sleep 1s
		
	fi
    
    
done

done


done



sed -i '/^$/d' ${r_log_jobs}PHASE_2_jobs.txt
while (( $(squeue -u adanguy | grep -f ${r_log_jobs}PHASE_2_jobs.txt | wc -l) >= 1 )) 
do    
    sleep 1s
done




rep_out=${r_log_PHASE}
rep_sortie1=${r_PHASE_output_format}
rep_sortie2=${r_intervals_mqs}
rep_sortie3=${r_PHASE_output_summary}
rep_sortie4=${r_graphs_PHASE}
rep_sortie5=${r_PHASE_time}

# Recherche de pb
beug=$(grep -Ril -e "segfault" -e "error" ${rep_out}*PHASE*)  
if [ ! -z "${beug}" ]
then
# Affichage (enlever extension)
beug2=$(echo ${beug} | sed "s|${rep_out}||g" | sed "s|.out||g" | sed "s|PHASE_||g")
echo ${beug2[@]} | tr ' ' '\n' | sort
echo ${beug2[@]} | tr ' ' '\n' | sort > ${r_log_PHASE}PHASE_fails.txt
# supression fichiers sorties (rq : le fichier log et sortie doivent s'appeler pareil sauf les repertoires et extensions)

for pb in ${beug2[@]} ;
do 
supression1=$(echo ${rep_sortie1}${pb}_PHASE_raw_outputs.txt)
supression2=$(echo ${rep_sortie2}${pb}_mqs.txt)
supression3=$(echo ${rep_sortie3}${pb}_PHASE_summary_outputs.txt)
supression4=$(echo ${rep_sortie4}${pb}.png)
supression5=$(echo ${rep_sortie5}${pb}_time.txt)
rm ${supression1}
rm ${supression2}
rm ${supression3}
rm ${supression4}
rm ${supression5}
done
fi




