#!/bin/bash

base=${1}



source ${base}





Rscript ${r_scripts}Data_embargo.R ${PLINK_initial}.map ${mks280k} ${r_amont}SNP_positions.txt ${r_PLINK}SNP_PLINK_4.txt ${r_published}SNP_positions_280k.txt

plink --file ${PLINK_initial} --recode --noweb --extract ${r_PLINK}SNP_PLINK_4.txt --out ${r_published}landraces_initial_genotyping_280k >>${r_log}Data_embargo.out


<<COMMENTS


repertoire=${2}

titre_resume=${r_PHASE}PHASE_summary_outputs.txt
titre_resume2=${r}common_SNP/PHASE/PHASE_summary_outputs.txt

titre_PHASE_summary_outputs2=${r_published}supplementary_file_S7_PHASE_summary_outputs.txt

titre_csre_genetic_map=${r_csre}csre_genetic_map.txt
titre_landraces_genetic_map=${r_maps_pop}landraces_genetic_maps.txt
titre_genetic_maps=${r_published}supplementary_file_S3_recombination_maps.txt

titre_hr=${r_HR}HR.txt
titre_hr2=${r_published}supplementary_file_S9_HRIs.txt
titre_cliques=${r_HR}cliques_HR.txt
titre_cliques2=${r_published}supplementary_file_S10_hot_windows.txt

titre_rho=${r_modele}chrregion_rho_random.txt
titre_rho2=${r}common_SNP/modele/chrregion_rho_random.txt

titre_lambda=${r_modele}chrregion_lambda_random.txt
titre_lambda2=${r}common_SNP/modele/chrregion_lambda_random.txt

titre_correlations_mixed_models=${r_published}supplementary_file_S8_correlations_mixed_models.txt

titre_landraces_4Mb=${r_mb}map_4mb_landraces.txt
titre_csre_4Mb=${r_csre}map_4mb_csre.txt
titre_cor_4Mb=${r_published}supplementary_file_S2_correlations_4Mb.txt
titre_maps_4Mb=${r_published}supplementary_file_S12_historical_and_meiotic_maps_4Mb.txt

titre_permutations_HR=${r_HR}permutations.txt
titre_overlap_HR=${r_HR}coloc_HR.txt
titre_overlap_HR_published=${r_published}supplementary_file_S11_overlap_HRIs.txt



titre_FST_SNP=${r_FST}chrregion_FST_SNP.txt
titre_FST_haplotypic_blocs=${r_FST}chrregion_FST_haplotypic_blocks.txt
titre_FST=${r_published}supplementary_file_S6_FST.txt

titre_SNP_positions=${r_amont}SNP_positions.txt
titre_chr_regions=${r_amont}Decoupage_chr_ble.tab
titre_correspondance_chr=${r_amont}Codes_chr.txt
titre_stab1=${r_published}supplementary_file_S5_region_limits.txt

titre_landraces=${r_landraces}landraces.txt
titre_stab2=${r_published}supplementary_file_S1_landraces_info.txt



##
v0=${repertoire}
v1=${titre_resume}
v2=${titre_resume2}
v3=${titre_PHASE_summary_outputs2}

v4=${titre_csre_genetic_map}
v5=${titre_landraces_genetic_map}
v6=${titre_genetic_maps}

v7=${titre_hr}
v8=${titre_hr2}
v9=${titre_cliques}
v10=${titre_cliques2}

v11=${titre_rho}
v12=${titre_rho2}

v13=${titre_lambda}
v14=${titre_lambda2}

v15=${titre_correlations_mixed_models}

v16=${titre_landraces_4Mb}
v17=${titre_csre_4Mb}
v18=${titre_cor_4Mb}
v19=${titre_maps_4Mb}

v20=${titre_permutations_HR}
v21=${titre_overlap_HR}


v22=${titre_overlap_HR_published}


v23=${titre_FST_SNP}
v24${titre_FST_haplotypic_blocs}
v25=${titre_FST}

v26=${titre_SNP_positions}
v27=${titre_chr_regions}
v28=${titre_correspondance_chr}
v29=${titre_stab1}

v30=${titre_landraces}
v31=${titre_stab2}



Rscript ${r_scripts}tabs_output.R ${v0} ${v1} ${v2} ${v3} ${v4} ${v5} ${v6} ${v7} ${v8} ${v9} ${v10} ${v11} ${v12} ${v13} ${v14} ${v15} ${v16} ${v17} ${v18} ${v19} ${v20} ${v21} ${v22} ${v23} ${v24} ${v25} ${v26} ${v27} ${v28} ${v29} ${v30} ${v31}


cp ${r_FST}pairwise_FST_matrix_haplotypic_blocks.txt ${r_published}pairwise_FST_matrix_haplotypic_blocks.txt
cp ${r_HR}relativ_lambda_HR_populations.txt ${r_published}supplementary_file_S4_intensity_under_HRIs.txt
cp ${r_FST}hapflk_tree_genome.txt ${r_published}HPAFLK_tree_genome_SNP.txt


titre_correlations_4mb_published=${titre_cor_4Mb}
titre_FST_published=${titre_FST}
titre_correlations_mixed_models_published=${titre_correlations_mixed_models}
titre_HR_published=${titre_hr2}
titre_cliques_published=${titre_cliques2}
titre_matrice_fst_blocs=${r_published}pairwise_FST_matrix_haplotypic_blocks.txt
titre_genetic_maps_published=${titre_genetic_maps}
titre_maps_4Mb=${titre_maps_4Mb}
titre_f2=${r_graphs}Figure2.tiff
titre_f3=${r_graphs}Figure3.tiff
titre_f4=${r_graphs}Figure4.tiff
titre_f5=${r_graphs}Figure5.tiff
titre_f6=${r_graphs}Figure6.tiff
titre_graph_significativity_boxplots_pairs_of_populations=${r_graphs}significativity_boxplots_pairs_of_populations.tiff
titre_slopes_SNP_specific_or_common=${r_graphs}slopes_SNP_specific_or_common.tiff
titre_overlap_HR_published=${r_published}overlap_HR_published.txt
titre_intensity_under_HR_published=${r_published}intensity_under_HR_published.txt
titre_correlations_mixed_models_published_common_SNP=${r}common_SNP/tabs/correlations_mixed_models_published.txt
titre_slopes_SNP_specific_or_common=${r_graphs}slopes_SNP_specific_or_common.tiff
titre_stab2=${r_published}supplementary_tab2.txt
titre_HAPFLK_tree_SNP=${r_published}HPAFLK_tree_genome_SNP.txt
titre_matrice_FST_blocs_png=${r_graphs}landraces_FST_matrix_haplotypes.png
titre_figure1=${r_graphs}Figure1.tiff


v1=${repertoire}
v2=${titre_correlations_4mb_published}
v3=${titre_FST_published}
v4=${titre_correlations_mixed_models_published}
v5=${titre_HR_published}
v6=${titre_cliques_published}
v7=${titre_matrice_fst_blocs}
v8=${titre_genetic_maps_published}
v9=${titre_maps_4Mb}
v10=${titre_overlap_HR_published}
v11=${titre_intensity_under_HR_published}
v12=${titre_f2}
v13=${titre_f3}
v14=${titre_f4}
v15=${titre_f5}
v16=${titre_f6}
v17=${titre_graph_significativity_boxplots_pairs_of_populations}
v18=${titre_correlations_mixed_models_published_common_SNP}
v19=${titre_slopes_SNP_specific_or_common}
v20=${titre_stab2}
v21=${titre_HAPFLK_tree_SNP}
v22=${titre_matrice_FST_blocs_png}
v23=${titre_figure1}



if [ "${repertoire}" != "common_SNP" ]
then

while (( $(squeue -u adanguy | grep -e "publish" | wc -l) >= 2 )) 
do    
    sleep 1s
done

fi



Rscript ${r_scripts}graphs_output.R ${v1} ${v2} ${v3} ${v4} ${v5} ${v6} ${v7} ${v8} ${v9} ${v10} ${v11} ${v12} ${v13} ${v14} ${v15} ${v16} ${v17} ${v18} ${v19} ${v20} ${v21} ${v22} ${v23}


COMMENTS
