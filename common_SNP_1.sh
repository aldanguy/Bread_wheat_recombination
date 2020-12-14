#!/bin/bash



# Estimate historical recombination rates on a common polymorphic SNP dataset



base=${1}
# base=/work/adanguy/these/pipeline/scripts/base_020820.sh

source ${base}


source ${base_common_SNP}

# rm -rf ${r}
source ${base_common_SNP}



# Define a common SNP dataset
job_out=${r_log}common_SNP_2.out
job_name=cSNP
# job10=$(sbatch -o ${job_out} -J ${job_name} --parsable ${r_scripts}common_SNP_2.sh ${base})

################################# Etape 5 : Historical recombination profiles of landraces


## SNP windows for PHASE
job_out=${r_log}windows_1_common_SNP.out
job_name=windows_1
# job11=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job10} --parsable ${r_scripts}windows_1.sh ${base_common_SNP})


## PHASE software and formatting of raw files
job_out=${r_log}PHASE_1_common_SNP.out
job_name=PHASE_1
# job12=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job11} --parsable ${r_scripts}PHASE_1.sh ${base_common_SNP}) # can take time, change for unlimitq


## Indentification of intervals for comparison of recombination profiles
job_out=${r_log}intervals_common_SNP.out
job_name=intervals
# job13=$(sbatch -o ${job_out} -J ${job_name} --mem=10G --dependency=afterok:${job12} --parsable ${r_scripts}intervals.sh ${base_common_SNP})


## Recombination rate within intervals
job_out=${r_log}intervals_PHASE_1_common_SNP.out
job_name=intervals_PHASE_1
# job14=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=10G --dependency=afterok:${job13} --parsable ${r_scripts}intervals_PHASE_1.sh ${base_common_SNP})


## Formatting and merging of data
job_out=${r_log}merge_common_SNP.out
job_name=merge
# job15=$(sbatch -o ${job_out} -J ${job_name} --mem=20G --dependency=afterok:${job14} --parsable ${r_scripts}merge.sh ${base_common_SNP})

################################# Etape 6 : data analysis

## Mixed model
job_name=modele
job_out=${r_log}modele.out
# job101=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=10G --dependency=afterok:${job15} --parsable ${r_scripts}modele.sh ${base_common_SNP})




## Highly recombining intervals (HR)
job_name=HR_1
job_out=${r_log}HR.out
# job16=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --dependency=afterok:${job15} --parsable ${r_scripts}HR_1.sh ${base_common_SNP})



## Heterogeneity of recombination landscape and Gini coefficients
job_name=gini
job_out=${r_log}gini.out
# job106=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --dependency=afterok:${job16} --parsable ${r_scripts}gini.sh ${base_common_SNP})




## Genetic maps averaged within 4 Mb windows
job_name=mb
size=4
job_out=${r_log}${size}mb.out
# job104=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --dependency=afterok:${job15} --parsable ${r_scripts}mb_1.sh ${base_common_SNP} ${size}) 
job104=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --parsable ${r_scripts}mb_1.sh ${base_common_SNP} ${size}) 


## Extraction of meiotic recombination maps from historical maps
job_name=specific_rec_rates
job_out=${r_log}population_specific_meiotic_rec_rates.out
# job107=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job15} --parsable ${r_scripts}population_specific_meiotic_rec_rates.sh ${base_common_SNP})



## Extraction of meiotic recombination maps from historical maps
job_name=published
job_out=${r_log}published.out
# job108=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --dependency=afterok:${job104}:${job16}:${job102}:${job101}:${job107} --parsable ${r_scripts}graphs_output.sh ${base_common_SNP} common_SNP)
job108=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=10G --dependency=afterok:${job104} --parsable ${r_scripts}graphs_output.sh ${base_common_SNP} common_SNP)

