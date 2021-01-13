#!/bin/bash



base=${1}

source ${base}
r_PLINK_base=${r_PLINK}
r_sources_base=${r_sources}
co_positions_base=${r_csre}csre_crossovers_positions.txt
pairwise_FST_matrix_haplotypic_blocks_base=${r_FST}pairwise_FST_matrix_haplotypic_blocks.txt 
chrregion_FST_haplotypic_blocks=${r_FST}chrregion_FST_haplotypic_blocks.txt
chrregion_FST_SNP=${r_FST}chrregion_FST_SNP.txt
hapflk_tree_genome=${r_FST}hapflk_tree_genome.txt



source ${base_common_SNP}
cp ${r_PLINK_base}*.map ${r_PLINK}
cp ${r_PLINK_base}*.ped ${r_PLINK}
cp ${co_positions_base} ${r_csre}
cp ${pairwise_FST_matrix_haplotypic_blocks_base} ${r_FST}pairwise_FST_matrix_haplotypic_blocks.txt 
cp ${chrregion_FST_haplotypic_blocks} ${r_FST}chrregion_FST_haplotypic_blocks.txt
cp ${chrregion_FST_SNP} ${r_FST}chrregion_FST_SNP.txt
cp ${hapflk_tree_genome} ${r_FST}hapflk_tree_genome.txt


Rscript ${r_scripts}common_SNP_2.R ${r_PLINK}WE.map ${r_PLINK}EE.map ${r_PLINK}WA.map ${r_PLINK}EA.map ${r_intervals}mqs_polym.txt






for p in ${pop[*]}
do

echo ${p}

plink --file ${r_PLINK}${p} --recode --extract ${r_intervals}mqs_polym.txt --noweb --out ${r_PLINK}${p} >>${r_log}common_SNP.out

rm ${r_PLINK}${p}.log
done

