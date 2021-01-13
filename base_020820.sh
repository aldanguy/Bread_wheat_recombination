#!/bin/bash



# General

date=020820

r_amont=/work/adanguy/these/pipeline/amont/

r_scripts=/work/adanguy/these/pipeline/scripts/

r=/work/adanguy/these/pipeline/${date}/


first=${r_scripts}first_020820.sh

base=${r_scripts}base_020820.sh

base_common_SNP=${r_scripts}base_${date}_common_SNP.sh

base_SNP_specific_vs_common=${r_scripts}base_${date}_SNP_specific_vs_common.sh


mkdir -p ${r}

cd ${r}



# Log


r_log=${r}log/
mkdir -p ${r_log}

r_log_jobs=${r}log/jobs/
mkdir -p ${r_log_jobs}

r_log_windows=${r}log/windows/
mkdir -p ${r_log_windows}

r_log_PHASE=${r}log/PHASE/
mkdir -p ${r_log_PHASE}

r_log_intervals=${r}log/intervals/
mkdir -p ${r_log_intervals}

r_log_csre=${r}log/csre/
mkdir -p ${r_log_csre}

r_log_HR=${r}log/HR/
mkdir -p ${r_log_HR}

r_log_mb=${r}log/mb/
mkdir -p ${r_log_mb}

r_log_alpha=${r}log/alpha/
mkdir -p ${r_log_alpha}

r_log_FST=${r}log/FST/
mkdir -p ${r_log_FST}

# Sources

r_sources=${r}sources/
mkdir -p ${r_sources}


# Work repositories


r_modele=${r}modele/
mkdir -p ${r_modele}

r_FST=${r}FST/
mkdir -p ${r_FST}

r_HR=${r}HR/
mkdir -p ${r_HR}

r_HR_permutations=${r}HR/permutations/
mkdir -p ${r_HR_permutations}

r_HR_alpha=${r}HR/alpha/
mkdir -p ${r_HR_alpha}

r_maps_pop=${r}maps_pops/
mkdir -p ${r_maps_pop}

r_mb=${r}mb/
mkdir -p ${r_mb}

r_gini=${r}gini/
mkdir -p ${r_gini}

r_PLINK=${r}PLINK/
mkdir -p ${r_PLINK}

r_PHASE=${r}PHASE/
mkdir -p ${r_PHASE}


r_PHASE_windows_description=${r}PHASE/windows/description/
mkdir -p ${r_PHASE_windows_description}

r_PHASE_windows_list=${r}PHASE/windows/list/
mkdir -p ${r_PHASE_windows_list}

r_PHASE_input_inp=${r}PHASE/input/inp/
mkdir -p ${r_PHASE_input_inp}

r_PHASE_input_recode_phase_inp=${r}PHASE/input/recode_phase_inp/
mkdir -p ${r_PHASE_input_recode_phase_inp}

r_PHASE_time=${r}PHASE/time/
mkdir -p ${r_PHASE_time}

r_PHASE_output_raw=${r}PHASE/output/raw/
mkdir -p ${r_PHASE_output_raw}

r_PHASE_output_format=${r}PHASE/output/format/
mkdir -p ${r_PHASE_output_format}

r_PHASE_output_summary=${r}PHASE/output/summary/
mkdir -p ${r_PHASE_output_summary}

r_intervals=${r}intervals/
mkdir -p ${r_intervals}

r_intervals_mqs=${r}intervals/mqs/
mkdir -p ${r_intervals_mqs}

r_intervals_PHASE=${r}intervals/PHASE/
mkdir -p ${r_intervals_PHASE}



r_csre=${r}csre/
mkdir -p ${r_csre}

# Graphes

r_graphs=${r}graphs/
mkdir -p ${r_graphs}

r_graphs_interpolation=${r}graphs/interpolation/
mkdir -p ${r_graphs_interpolation}


r_graphs_PHASE=${r}graphs/PHASE/
mkdir -p ${r_graphs_PHASE}

r_published=${r}tabs/
mkdir -p ${r_published}

r_scaffolds=${r}scaffolds/
mkdir -p ${r_scaffolds}

r_landraces=${r}landraces/
mkdir -p ${r_landraces}

r_annotation=${r}annotation/
mkdir -p ${r_annotation}


##### Fichiers amont


scaffolds_input=${r_amont}161010_Chinese_Spring_v1.0_pseudomolecules_AGP.tsv

accessions=${r_amont}accessions_4403_alice.csv

structure=${r_amont}structure_tot_def.txt

K4=${r_amont}landrace8741hap_k4_rep5_f_f-ord-ord-ord

K8=${r_amont}landrace8741hap_k8_rep2_f_f-ord-ord-ord

dis=${r_amont}landrace632_8741hap.dis

chr_names=${r_amont}Codes_chr.txt

chr_tab=${r_amont}Decoupage_chr_ble.tab

PLINK_initial=${r_amont}BW_261K_wp3

haplotypic_blocks=${r_amont}landrace632_8741hap_sophie.var

haplotypic_blocks_limits=${r_amont}stat_haplo8741.txt

SNP_positions=${r_amont}SNP_positions.txt

amont_annotation=${r_amont}iwgsc_refseqv1.0_HighConf_UTR_2017May05.gff3

csre_genotyping=${r_amont}csre_genotyping.txt


mks280k=${r_amont}Rimbert_2018_tabs2.txt


## Variables

pop=(WE EE WA EA)

chr=(1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 7B 7D)

features=(gene five_prime_UTR exon three_prime_UTR)

nb_jobs_allowed=300

X=10

x=10

nb_posterior=$((X * 100))


rate_het_max=0.05
rate_NA_max=0.1
rate_maf_min=0.03


prior=${r_amont}prior_default.txt

burn_in=100


subset_HR=1000

nb_perm_HR=1000

module purge
# module load system/R-3.4.3
module load system/R-3.6.2
module load system/pandoc-2.1.3
module load bioinfo/bedtools-2.27.1
module load bioinfo/bedops-v2.4.35
module load bioinfo/tabix-0.2.5
module load system/Python-3.7.4

