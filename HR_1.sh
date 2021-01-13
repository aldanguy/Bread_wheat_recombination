#!/bin/bash

base=${1}


source ${base}


# With marges

borders=NA

job_name=HR_border

job_out=${r_log}HR_borders.out



job=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --parsable ${r_scripts}HR_2.sh ${base} ${borders})



while (( $(squeue -u adanguy | grep ${job} | wc -l) >= 1 )) 
do    
    sleep 1s
done




cp ${r_HR}permutations.txt ${r_HR}permutations_borders.txt

cp ${r_HR}coloc_HR.txt ${r_HR}coloc_HR_borders.txt 

cp ${r_HR}cliques_HR_${subset_HR}.txt ${r_HR}cliques_HR_${subset_HR}_borders.txt

cp ${r_HR}cliques_HR.txt ${r_HR}cliques_HR_borders.txt

cp ${r_graphs}coloc_HR.png ${r_graphs}coloc_HR_borders.png 

cp ${r_graphs}coloc_HR_distances.png ${r_graphs}coloc_HR_distances_borders.png

cp ${r_graphs}relative_lambda_HR_populations_3.png ${r_graphs}relative_lambda_HR_populations_3_borders.png 

cp ${r_graphs}Figure5.tiff ${r_graphs}Figure5_borders.tiff

cp ${r_HR}relativ_lambda_HR_populations.txt ${r_HR}relativ_lambda_HR_populations_borders.txt


## Without marges
borders=0

job_name=HR_no_border

job_out=${r_log}HR.out



job=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --parsable ${r_scripts}HR_2.sh ${base} ${borders})

while (( $(squeue -u adanguy | grep ${job} | wc -l) >= 1 )) 
do    
    sleep 1s
done

