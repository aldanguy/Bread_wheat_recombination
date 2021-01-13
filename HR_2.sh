#!/bin/bash

base=${1}


source ${base}

borders=${2}


Rscript ${r_scripts}HR_detection.R ${r_PHASE}PHASE_summary_outputs.txt ${chr_names} ${r_HR}HR.txt ${r_graphs}distribution_lambda.png ${r_graphs}ratio_HR_size.png ${r_graphs}relationship_lambda_size.png ${r_graphs}physical_size_HR.png
 

 
rm ${r_HR_permutations}*
rm ${r_log_HR}*



n=0

while (( ${n} <= ${nb_perm_HR} ))
do
	
	nb_jobs_squeue=$(squeue -u adanguy | wc -l)

	
	
	if (( ${nb_jobs_squeue} <= ${nb_jobs_allowed} )) 
	
	then
    
     echo ${n}
  
	

    job_out=${r_log_HR}permutation_${n}.out
    job_name=${n}_p
	    
    job_perm=$(sbatch -o ${job_out} -J ${job_name} --mem=40G --parsable ${r_scripts}HR_permutation.sh ${base} ${borders} ${n} )
    
    
    echo "${job_perm}" >> ${r_log_jobs}permutation_jobs.txt

		
    n=$((n +1))
		
	else
	
		sleep 1s
		
	fi
    
    
done

sed -i '/^$/d' ${r_log_jobs}permutation_jobs.txt
while (( $(squeue -u adanguy | grep -f ${r_log_jobs}permutation_jobs.txt | wc -l) >= 1 )) 
do    
    sleep 1s
done

beug=$(grep -Ril -e "segfault" -e "error" ${r_log_HR})
if [ ! -z "${beug}" ]
then
    # Affichage (enlever extension)
    beug2=$(echo ${beug} | sed "s|${r_log_HR}||g" | sed "s|.out||g")
    echo ${beug2[@]} | tr ' ' '\n' | sort
    echo ${beug2[@]} | tr ' ' '\n' | sort > ${r_log_HR}HR_fails.txt
    for pb in ${beug2[@]} ;
        do 
        supression1=$(echo ${r_HR_permutations}${pb}.txt)
        rm ${supression1}
    done
fi






k=0

for f in ${r_HR_permutations}permutation_*.txt
do 

    if ((k==0 ))
        
    then 
        
        cat ${f} > ${r_HR}permutations.txt
            # rm ${f}
        
    else
        
        tail -n+2 ${f} >> ${r_HR}permutations.txt
            # rm ${f}

    fi
        
    k=$((k +1))
        
          
done

Rscript ${r_scripts}HR_detection_${subset_HR}.R ${r_HR}HR.txt ${r_HR}HR_${subset_HR}.txt 

Rscript ${r_scripts}HR_coloc.R reference ${r_HR}HR.txt ${r_PHASE}PHASE_summary_outputs.txt ${r_HR}coloc_HR.txt ${borders} ${r_HR}cliques_HR.txt

Rscript ${r_scripts}HR_coloc.R ${subset_HR} ${r_HR}HR_${subset_HR}.txt ${r_PHASE}PHASE_summary_outputs.txt ${r_HR}coloc_HR_${subset_HR}.txt ${borders} ${r_HR}cliques_HR_${subset_HR}.txt

COMMENTS

Rscript ${r_scripts}HR_graphs.R ${r_HR}permutations.txt ${r_HR}cliques_HR.txt ${r_HR}HR.txt ${r_HR}coloc_HR.txt ${r_PHASE}PHASE_summary_outputs.txt ${r_HR}cliques_HR_${subset_HR}.txt ${r_graphs}coloc_HR.png ${r_graphs}coloc_HR_distances.png ${r_graphs}Figure5.tiff ${r_graphs}relative_lambda_HR_populations_2.png ${r_HR}relativ_lambda_HR_populations.txt ${r_HR}relativ_lambda_HR_combination_populations.txt ${r_graphs}relative_lambda_HR_populations_3.png



