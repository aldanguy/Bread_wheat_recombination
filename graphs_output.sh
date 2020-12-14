#!/bin/bash

base=${1}


source ${base}





repertoire=${2}

titre_resume=${r_PHASE}PHASE_summary_outputs.txt
titre_PHASE_summary_outputs2=${r_published}PHASE_summary_outputs_published.txt

titre_csre_genetic_map=${r_sources}csre_genetic_map.txt
titre_landraces_genetic_map=${r_maps_pop}landraces_genetic_maps.txt
titre_genetic_maps=${r_published}genetic_maps_published.txt

titre_hr=${r_HR}HR.txt
titre_hr2=${r_published}HR_published.txt
titre_cliques=${r_HR}cliques_HR.txt
titre_cliques2=${r_published}cliques_published.txt

titre_rho=${r_modele}chrregion_rho_random.txt
titre_lambda=${r_modele}chrregion_lambda_random.txt
titre_correlations_mixed_models=${r_published}correlations_mixed_models_published.txt

titre_landraces_4Mb=${r_mb}map_4mb_landraces.txt
titre_csre_4Mb=${r_csre}map_4mb_csre.txt
titre_cor_4Mb=${r_published}correlations_4mb_published.txt
titre_maps_4Mb=${r_published}historical_and_meiotic_maps_4Mb_published.txt



titre_csre_genotyping=${r_sources}csre_genotyping.txt
titre_csre_genotyping2=${r_published}csre_genotyping_matrix_published.txt


titre_FST_SNP=${r_FST}chrregion_FST_SNP.txt
titre_FST_haplotypic_blocs=${r_FST}chrregion_FST_haplotypic_blocks.txt
titre_FST=${r_published}FST_published.txt

titre_SNP_positions=${r_sources}SNP_positions.txt
titre_chr_regions=${r_amont}Decoupage_chr_ble.tab
titre_correspondance_chr=${r_amont}Codes_chr.txt
titre_stab1=${r_published}supplementary_tab1.txt

titre_landraces=${r_sources}landraces.txt
titre_stab2=${r_published}supplementary_tab2.txt



##
v0=${repertoire}
v1=${titre_resume}
v2=${titre_PHASE_summary_outputs2}

v3=${titre_csre_genetic_map}
v4=${titre_landraces_genetic_map}
v5=${titre_genetic_maps}

v6=${titre_hr}
v7=${titre_hr2}
v8=${titre_cliques}
v9=${titre_cliques2}

v10=${titre_rho}
v11=${titre_lambda}
v12=${titre_correlations_mixed_models}

v13=${titre_landraces_4Mb}
v14=${titre_csre_4Mb}
v15=${titre_cor_4Mb}
v16=${titre_maps_4Mb}


v17=${titre_csre_genotyping}
v18=${titre_csre_genotyping2}


v19=${titre_FST_SNP}
v20=${titre_FST_haplotypic_blocs}
v21=${titre_FST}

v22=${titre_SNP_positions}
v23=${titre_chr_regions}
v24=${titre_correspondance_chr}
v25=${titre_stab1}

v26=${titre_landraces}
v27=${titre_stab2}



Rscript ${r_scripts}tabs_output.R ${v0} ${v1} ${v2} ${v3} ${v4} ${v5} ${v6} ${v7} ${v8} ${v9} ${v10} ${v11} ${v12} ${v13} ${v14} ${v15} ${v16} ${v17} ${v18} ${v19} ${v20} ${v21} ${v22} ${v23} ${v24} ${v25} ${v26} ${v27}


cp ${r_FST}pairwise_FST_matrix_haplotypic_blocks.txt ${r_published}pairwise_FST_matrix_haplotypic_blocks.txt



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
titre_f6=${r_graphs}Figure6.tiff
titre_graph_significativity_boxplots_pairs_of_populations=${r_graphs}significativity_boxplots_pairs_of_populations.tiff

v1=${repertoire}
v2=${titre_correlations_4mb_published}
v3=${titre_FST_published}
v4=${titre_correlations_mixed_models_published}
v5=${titre_HR_published}
v6=${titre_cliques_published}
v7=${titre_matrice_fst_blocs}
v8=${titre_genetic_maps_published}
v9=${titre_maps_4Mb}
v10=${titre_f2}
v11=${titre_f3}
v12=${titre_f4}
v13=${titre_f6}
v14=${titre_graph_significativity_boxplots_pairs_of_populations}


Rscript ${r_scripts}graphs_output.R ${v1} ${v2} ${v3} ${v4} ${v5} ${v6} ${v7} ${v8} ${v9} ${v10} ${v11} ${v12} ${v13} ${v14}
