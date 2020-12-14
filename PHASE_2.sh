#!/bin/bash

SECONDS=0


base=${1}

source ${base}

p=${2}

c=${3}

f=${4}

ID=${p}_${c}_${f}

recom=${r_PHASE_output_raw}${ID}_recom
monitor=${r_PHASE_output_raw}${ID}_monitor
inp=${r_PHASE_input_inp}${ID}.inp
out=${r_PHASE_output_format}${ID}
window=${r_PHASE_windows_description}w_${p}_${c}.txt




###### Conversion au format fastPHASE

# plink --file ${r_PLINK}${p} --recode-fastphase --extract ${r_PHASE_windows_list}w_${ID}.txt --noweb --out ${r_PHASE_input_recode_phase_inp}${ID} >>${r_log_PHASE}PHASE_${p}_${c}_${f}.out

# rm ${r_PHASE_input_recode_phase_inp}${ID}.log

###### Conversion pour format PHASE


# Nombre de marqueurs 
N=$(head -n2 ${r_PHASE_input_recode_phase_inp}${ID}.recode.phase.inp | tail -n1)





# Variable contenant autant de "S" que de marqueurs conserves apres filtres de PLINK
myvar=`perl -e "print 'S' x ${N};"`

sed '4i\'${myvar}'\' ${r_PHASE_input_recode_phase_inp}${ID}.recode.phase.inp | sed 's/ ID //g' > ${r_PHASE_input_inp}${ID}.inp





##### Etape 3 PHASE




# /save/servin/bin/PHASE -S${f} -X${X} -r${prior} -k999 -x${x} ${r_PHASE_input_inp}${ID}.inp ${r_PHASE_output_raw}${ID} 100 1 ${burn_in} >>${r_log_PHASE}PHASE_${p}_${c}_${f}.out 2>${r_log_PHASE}PHASE_${p}_${c}_${f}.err
# if outputs from PHASE show no variation or problem in PAC likelihood, change seed=seed+1 and add 100 to burn_in -> 100 1 200


# Mise en forme




# Nb intervalles
ni=$((${N} -1))


nb_underflows=$(cat ${r_log_PHASE}PHASE_${p}_${c}_${f}.err | grep -e "underflow problem in computation of logFDLSProb" | wc -l) ;
    
nb_fails=$(cat ${r_log_PHASE}PHASE_${p}_${c}_${f}.err | grep -e "failed to find decent estimate of recombination parameters" | wc -l) ;




cat ${recom} |
cut -f1 -d" "|
tail -n+2 |
while read line;
do   
    for i in $(seq 1 ${ni}) ;
    do     echo "${line}";
    done;
done > ${out}_rho.txt


cat ${recom} | cut -d" " -f2- | tail -n+2 | awk '{for (i=1; i<=NF; i++) a[i,NR]=$i
        max=(max<NF?NF:max)}
        END {for (i=1; i<=max; i++)
              {for (j=1; j<=NR; j++) 
                  printf "%s%s", a[i,j], (j==NR?RS:FS)
              }
        }' | tr ' ' '\n' > ${out}_lambda.txt



cat ${inp} |
tail -n +3 |
head -n1 |
cut -f2- -d" " |
tr -s ' '  '\n' |
awk -F $' ' ' BEGIN {OFS = FS} NR > 1 { print prev, $1 } { prev = $0 } END { print prev }' |
head -n -1 |
while read line; 
do   
    for i in $(seq 1 ${nb_posterior}) ; 
    do     echo "${line}";   
    done; 
done  > ${out}_pos.txt



paste -d" " ${out}_pos.txt ${out}_rho.txt ${out}_lambda.txt |
tr ' ' '\t' |
sed "s/^/${ID}\t/" |
sed "s/^/${c}\t/" |
sed "s/^/${p}\t/" |
sed "s/$/\t${nb_underflows}/" |
sed "s/$/\t${nb_fails}/" > ${out}_PHASE_raw_outputs_temp.txt


echo -e "population\tchr\tID\tposSNPlpop\tposSNPrpop\trho\tlambda\tnb_underflows\tnb_fails" > ${out}_PHASE_raw_outputs.txt

cat ${out}_PHASE_raw_outputs_temp.txt >> ${out}_PHASE_raw_outputs.txt


rm ${out}_pos.txt
rm ${out}_rho.txt
rm ${out}_lambda.txt
rm ${out}_PHASE_raw_outputs_temp.txt



### Ajout d'infos



Rscript ${r_scripts}PHASE_2.R ${monitor} ${window} ${r_sources}scaffolds.txt ${out}_PHASE_raw_outputs.txt ${r_sources}SNP_positions.txt ${out}_PHASE_raw_outputs.txt ${r_intervals_mqs}${ID}_mqs.txt ${r_PHASE_output_summary}${ID}_PHASE_summary_outputs.txt ${r_graphs_PHASE}${ID}.png


# Info sur la fenetre

PAC_slope=$(head ${out}_PHASE_raw_outputs.txt | cut -f16 | tail -n1)



# echo -e "${ID}\t${SECONDS}\t${nb_underflows}\t${nb_fails}\t${PAC_slope}" > ${r_PHASE_time}${ID}_time.txt


