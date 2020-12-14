#!/bin/bash

## Figure 1 of article


base=${1}


source ${base}

Rscript ${r_scripts}landraces_graphs.R ${r_sources}landraces.txt ${r_FST}hapflk_tree_genome.txt ${r_FST}pairwise_FST_matrix_haplotypic_blocks.txt ${dis} ${r_graphs}landraces1.png ${r_graphs}Figure1.tiff ${r_graphs}landraces_FST_matrix_haplotypes.png ${r_graphs}landraces_hclust4.png ${r_sources}tabs2.txt

