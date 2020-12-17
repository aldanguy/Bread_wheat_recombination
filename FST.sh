#!/bin/bash


## Compute Fst either with haplotypic blocks or SNP data


base=${1}


source ${base}



rm ${r_FST}*

############# With haplotypic blocks



module purge

module load system/R-3.4.3_bis



Rscript ${r_scripts}FST_blocks.R ${chr_tab} ${chr_names} ${haplotypic_blocks} ${haplotypic_blocks_limits} ${r_landraces}landraces.txt ${r_FST}chrregion_FST_haplotypic_blocks.txt ${r_FST}pairwise_FST_matrix_haplotypic_blocks.txt



############# With SNP 
source ${base}

# Pairwise Reynolds distance

source /save/servin/Envs/hapflkdev/bin/activate

hapflk --bfile ${r_landraces}landraces_final >>${r_log}FST.out

deactivate


cp ${r}hapflk_reynolds.txt ${r_FST}hapflk_reynolds_genome.txt
cp ${r}hapflk_tree.txt ${r_FST}hapflk_tree_genome.txt

for f in ${r_PLINK}*chrregion*; 
do

chrregion=$(echo ${f} | sed "s|${r_PLINK}||g" | sed "s|.txt||g" | sed "s|chrregion_||g")

echo ${chrregion}


plink --file ${r_landraces}landraces_final --recode --noweb --extract ${f} --out ${r_FST}landraces_${chrregion} >>${r_log}FST.out

plink --file ${r_FST}landraces_${chrregion} --noweb --make-bed --out ${r_FST}landraces_${chrregion} >>${r_log}FST.out

source /save/servin/Envs/hapflkdev/bin/activate

hapflk --bfile ${r_FST}landraces_${chrregion} >>${r_log}FST.out

deactivate





cp ${r}hapflk_reynolds.txt ${r_FST}hapflk_reynolds_${chrregion}.txt


Rscript ${r_scripts}FST_SNP.R ${chrregion} ${r_FST}hapflk_reynolds_${chrregion}.txt ${f} ${r_FST}chrregion_FST_SNP.txt


rm ${r}hapflk*
rm ${r_FST}landraces_${chrregion}*


done

