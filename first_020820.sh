#!/bin/bash
#SBATCH -o /work/adanguy/these/pipeline/scripts/pipeline.out
#SBATCH -J first
#SBATCH --time=00:10:00


################################# Etape 1 : repositories


base=/work/adanguy/these/pipeline/scripts/base_020820.sh

source ${base}

# rm -rf ${r}

source ${base}



################################# Etape 2 : sources files


## scaffolds limits
job_out=${r_log}scaffolds.out
job_name=scaffolds
# job1=$(sbatch -o ${job_out} -J ${job_name} --time=00:30:00 --parsable ${r_scripts}scaffolds.sh ${base})


## Passeport data about landraces
job_out=${r_log}landraces.out
job_name=landraces
# job2=$(sbatch -o ${job_out} -J ${job_name} --time=00:30:00 --parsable ${r_scripts}landraces.sh ${base})


## Positions of SNP
job_out=${r_log}SNP_positions.out
job_name=SNP_positions
# job3=$(sbatch -o ${job_out} -J ${job_name} --time=00:30:00 --parsable ${r_scripts}SNP_positions.sh ${base})
# This script is skipped in public pipeline because of management of not public data. Not necessary to run further scripts.

## CsRe genotyping, crossovers counting and empirical genetic map 
job_out=${r_log_csre}csre.out
job_name=csre
# job4=$(sbatch -o ${job_out} -J ${job_name} --time=48:00:00 --mem=50G --parsable ${r_scripts}csre.sh ${base})
# Within this bash script, one R script is skipped in public pipeline because of management of not public data. This skipped script is not necessary to run further scripts.


################################# Etape 3 : CsRe Bayesian recombination profile


## CsRe Bayesian genetic map
job_out=${r_log_csre}csre_bayesian.out
job_name=csre_bayesian
# job5=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --dependency=afterok:${job4} --parsable ${r_scripts}csre_bayesian.sh ${base} )


## Interpolation of genetic distances of unmapped SNP according CsRe Bayesian map
job_out=${r_log}interpolation.out
job_name=interpolation
# job6=$(sbatch -o ${job_out} -J ${job_name} --time=00:30:00 --dependency=afterok:${job5} --parsable ${r_scripts}interpolation.sh ${base})


## Update of PLINK files (physical and genetic positions of SNP)
job_out=${r_log}PLINK.out
job_name=PLINK
# job7=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --dependency=afterok:${job6} --parsable ${r_scripts}PLINK.sh ${base})


################################# Etape 4 : Identification of populations of landraces

## 4 diverging populations of landraces
job_out=${r_log}landraces_structuration.out
job_name=landraces_structuration
# job8=$(sbatch -o ${job_out} -J ${job_name} --time=00:30:00 --dependency=afterok:${job7} --parsable ${r_scripts}landraces_structuration.sh ${base})


## Filtering on Minor Allelic Frequency within each population
job_out=${r_log}MAF_filtering.out
job_name=MAF
# job9=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --dependency=afterok:${job8} --parsable ${r_scripts}MAF_filtering.sh ${base} )

## FST of landraces populations
job_out=${r_log}FST.out
job_name=FST
# job10=$(sbatch -o ${job_out} -J ${job_name} --mem=10G --dependency=afterok:${job9} --parsable ${r_scripts}FST.sh ${base})

## Representation for landraces
job_name=graphe
job_out=${r_log}landraces_graphe.out
# job100=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job10} --parsable ${r_scripts}landraces_graphs.sh ${base})



################################# Etape 5 : Historical recombination profiles of landraces


## SNP windows for PHASE
job_out=${r_log}windows_1.out
job_name=windows_1
# job11=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job9} --parsable ${r_scripts}windows_1.sh ${base})


## PHASE software and formatting of raw files
job_out=${r_log}PHASE_1.out
job_name=PHASE_1
# job12=$(sbatch -o ${job_out} -J ${job_name} -p unlimitq --dependency=afterok:${job11} --parsable ${r_scripts}PHASE_1.sh ${base})


## Indentification of intervals for comparison of recombination profiles
job_out=${r_log}intervals.out
job_name=intervals
# job13=$(sbatch -o ${job_out} -J ${job_name} --mem=10G --dependency=afterok:${job12} --parsable ${r_scripts}intervals.sh ${base})



## Recombination rate within intervals
job_out=${r_log}intervals_PHASE_1.out
job_name=intervals_PHASE_1
# job14=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=10G --dependency=afterok:${job13} --parsable ${r_scripts}intervals_PHASE_1.sh ${base})


## Formatting and merging data
job_out=${r_log}merge.out
job_name=merge
# job15=$(sbatch -o ${job_out} -J ${job_name} --mem=20G --dependency=afterok:${job14} --parsable ${r_scripts}merge.sh ${base})

################################# Etape 6 : data analysis


## Mixed model
job_name=modele
job_out=${r_log}modele.out
# job101=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=10G --dependency=afterok:${job15} --parsable ${r_scripts}modele.sh ${base})




## Check that model results does not depend on polymorphic SNP 
job_name=cSNP
job_out=${r_log}common_SNP.out
# job102=$(sbatch -o ${job_out} -J ${job_name} -p workq --time=00:10:00 --dependency=afterok:${job10} --parsable ${r_scripts}common_SNP_1.sh ${base})
job102=$(sbatch -o ${job_out} -J ${job_name} -p workq --time=00:10:00 --parsable ${r_scripts}common_SNP_1.sh ${base})


## Highly recombining intervals (HR)
job_name=HR_1
job_out=${r_log}HR.out
# job16=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --dependency=afterok:${job15} --parsable ${r_scripts}HR_1.sh ${base})
job16=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --parsable ${r_scripts}HR_1.sh ${base})


## Position of gene features and HR
job_name=at
job_out=${r_log}annotation.out
# job103=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job16} --parsable ${r_scripts}annotations.sh ${base})







## Genetic maps averaged within 4 Mb windows
job_name=mb
size=4
job_out=${r_log}${size}mb.out
# job104=$(sbatch -o ${job_out} -J ${job_name} --mem=50G --dependency=afterok:${job15} --parsable ${r_scripts}mb_1.sh ${base} ${size}) 

## Relation between CsRe recombination and density of landraces HR
job_name=d_i
job_out=${r_log}csre_and_HR.out
# job105=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job16}:${job104} --parsable ${r_scripts}csre_and_HR.sh ${base})





## Heterogeneity of recombination landscape and Gini coefficients
job_name=gini
job_out=${r_log}gini.out


## Extraction of meiotic recombination maps from historical maps
job_name=specific_rec_rates
job_out=${r_log}population_specific_meiotic_rec_rates.out
# job107=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job15} --parsable ${r_scripts}population_specific_meiotic_rec_rates.sh ${base})


## Extraction of meiotic recombination maps from historical maps
job_name=published
job_out=${r_log}published.out
# job108=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --dependency=afterok:${job104}:${job102}:${job16}:${job102}:${job101}:${job107} --parsable ${r_scripts}graphs_output.sh ${base} all)
job108=$(sbatch -o ${job_out} -J ${job_name} -p workq --mem=50G --parsable ${r_scripts}published.sh ${base} all)


## comparison of recombination rates with different SNP dataset
job_name=SNP_specific_vs_common
job_out=${r_log}SNP_specific_vs_common.out
# job109=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job12}:${job102} --parsable ${r_scripts}SNP_specific_vs_common_1.sh ${base})


## comparison of recombination rates with different SNP dataset
job_name=LD_window
job_out=${r_log}LD_window.out
# job110=$(sbatch -o ${job_out} -J ${job_name} -p workq --dependency=afterok:${job15} --parsable ${r_scripts}select_window_for_DL_graph.sh ${base})







