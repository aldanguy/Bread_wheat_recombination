#!/bin/bash

base=${1}


source ${base}

size=${2}

step=${2}



size2=$(awk -v size="$size" 'BEGIN {print (size*1000000)}')

step2=$(awk -v step="$step" 'BEGIN {print (step*1000000)}')




titre=map_${size}mb


map=${r_mb}${titre}.txt


rm ${r_mb}*${titre}*



Rscript ${r_scripts}mb.R ${step2} ${size2} ${map}





for c in ${chr[*]}
do


for p in ${pop[*]}
do



    job_out=${r_log_mb}${titre}_landraces_${p}_${c}.out

    job_name=${p}_${c}_mb
    
    job=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --parsable ${r_scripts}mb_2.sh ${base} ${p} ${c} ${map} ${titre}_landraces_${p}_${c})
    
    echo "${job}" >> ${r_log_jobs}map_${size}mb_2_jobs.txt
    


done
done






Rscript ${r_scripts}mb_csre.R ${r_csre}csre_crossovers_positions.txt ${r_csre}csre_genetic_map.txt ${chr_tab} ${nb_posterior} ${map} ${r_csre}${titre}_csre.txt


sed -i '/^$/d' ${r_log_jobs}map_${size}mb_2_jobs.txt
while (( $(squeue -u adanguy | grep -f ${r_log_jobs}map_${size}mb_2_jobs.txt | wc -l) >= 1 )) 
do    
    sleep 1s
done




rep_out=${r_log_mb}
rep_sortie=${r_mb}
# Recherche de pb
beug=$(grep -Ril -e "segfault" -e "error" ${rep_out}*${titre}*)  
if [ ! -z "${beug}" ]
then
beug2=$(echo ${beug} | sed "s|${rep_out}||g" | sed "s|.out||g")
# supression fichiers sorties (rq : le fichier log et sortie doivent s'appeler pareil sauf les repertoires et extensions)
supression=$(echo ${beug} | sed "s|${rep_out}|${rep_sortie}|g" | sed "s|.out|.txt|g")
echo ${beug2}
rm ${supression}
fi





k=0

for f in ${r_mb}${titre}_landraces_*.txt
do 

    if ((k==0 ))
        
    then 
        
        cat ${f} > ${r_mb}${titre}_landraces.txt
            rm ${f}
        
    else
        
        tail -n+2 ${f} >> ${r_mb}${titre}_landraces.txt
            rm ${f}

    fi
        
    k=$((k +1))
        
          
done

rm ${map}
