

#### Graphics

Sys.time()
cat("\n\nGraphs_output.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)



suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(qgraph))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(gridGraphics))
suppressPackageStartupMessages(library(nlme))


variables <- commandArgs(trailingOnly=TRUE)

cat("\n Graphics\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

repertoire <- variables[1]
titre_correlations_4mb_published <- variables[2]
titre_FST_published <- variables[3]
titre_correlations_mixed_models_published <- variables[4]
titre_HR_published <- variables[5]
titre_cliques_published<- variables[6]
titre_matrice_FST_blocs<- variables[7]
titre_genetic_maps_published<- variables[8]
titre_maps_4Mb <- variables[9]
titre_overlap_HR_published <- variables[10]
titre_intensity_under_HR_published <- variables[11]
titre_f2 <- variables[12]
titre_f3 <- variables[13]
titre_f4 <- variables[14]
titre_f5 <- variables[15]
titre_f6 <- variables[16]
titre_graph_significativity_boxplot_pairs_of_populations <- variables[17]
titre_correlations_mixed_models_published_common_SNP <- variables[18]
titre_slopes_SNP_specific_or_common<- variables[19]
titre_landraces_published <- variables[20]
titre_HAPFLK_tree_SNP <- variables[21]
titre_matrice_FST_blocs_png <-  variables[22]
titre_figure1 <-  variables[23]





# 
# repertoire <- "all"                                                                                          
# titre_correlations_4mb_published <-  "/work/adanguy/these/pipeline/020820/tabs/correlations_4mb_published.txt"                      
# titre_FST_published <-"/work/adanguy/these/pipeline/020820/tabs/FST_published.txt"                                   
# titre_correlations_mixed_models_published <-"/work/adanguy/these/pipeline/020820/tabs/correlations_mixed_models_published.txt"             
# titre_HR_published <- "/work/adanguy/these/pipeline/020820/tabs/HR_published.txt"                                    
# titre_cliques_published<- "/work/adanguy/these/pipeline/020820/tabs/cliques_published.txt"                               
# titre_matrice_FST_blocs<- "/work/adanguy/these/pipeline/020820/tabs/pairwise_FST_matrix_haplotypic_blocks.txt"           
# titre_genetic_maps_published<-  "/work/adanguy/these/pipeline/020820/tabs/genetic_maps_published.txt"                          
# titre_maps_4Mb <- "/work/adanguy/these/pipeline/020820/tabs/historical_and_meiotic_maps_4Mb_published.txt"       
# titre_overlap_HR_published <- "/work/adanguy/these/pipeline/020820/tabs/overlap_HR_published.txt"                            
# titre_intensity_under_HR_published <-  "/work/adanguy/these/pipeline/020820/tabs/intensity_under_HR_published.txt"                    
# titre_f2 <- "/work/adanguy/these/pipeline/020820/graphs/Figure2.tiff"                                      
# titre_f3 <- "/work/adanguy/these/pipeline/020820/graphs/Figure3.tiff"                                      
# titre_f4 <-"/work/adanguy/these/pipeline/020820/graphs/Figure4.tiff"                                      
# titre_f5 <- "/work/adanguy/these/pipeline/020820/graphs/Figure5.tiff"                                      
# titre_f6 <- "/work/adanguy/these/pipeline/020820/graphs/Figure6.tiff"                                      
# titre_graph_significativity_boxplot_pairs_of_populations <-  "/work/adanguy/these/pipeline/020820/graphs/significativity_boxplots_pairs_of_populations.tiff"
# titre_correlations_mixed_models_published_common_SNP <- "/work/adanguy/these/pipeline/020820/common_SNP/tabs/correlations_mixed_models_published.txt"  
# titre_slopes_SNP_specific_or_common<- "/work/adanguy/these/pipeline/020820/graphs/slopes_SNP_specific_or_common.tiff"                
# titre_landraces_published <- "/work/adanguy/these/pipeline/020820/tabs/supplementary_tab2.txt"                              
# titre_HAPFLK_tree_SNP <- "/work/adanguy/these/pipeline/020820/tabs/HPAFLK_tree_genome_SNP.txt"                          
# titre_matrice_FST_blocs_png <- "/work/adanguy/these/pipeline/020820/graphs/landraces_FST_matrix_haplotypes.png"               
# titre_figure1 <-  "/work/adanguy/these/pipeline/020820/graphs/Figure1.tiff" 
# 
# 







   titre_correlations_4mb_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/correlations_4mb_published.txt"
   titre_FST_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/FST_published.txt"
   titre_correlations_mixed_models_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/supplementary_file_S8_correlations_mixed_models.txt"
   titre_HR_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/HR_published.txt"
   titre_cliques_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/titre_cliques_published.txt"
   titre_matrice_FST_blocs <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/pairwise_FST_matrix_haplotypic_blocks.txt"
   titre_genetic_maps_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/genetic_maps_published.txt"
   titre_maps_4Mb <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/historical_and_meiotic_maps_4Mb_published.txt"
   titre_f3 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure3.tiff"
   titre_f2 <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure2.tiff"
   titre_f4 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure4.tiff"
   titre_f6 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure6.tiff"
   titre_f6 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes_SNP_common/Figure6.tiff"
   titre_slopes_SNP_specific_or_common="/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/slopes_SNP_specific_or_common.tiff"
   titre_f5 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes_SNP_common/Figure5.tiff"
   
  
   titre_overlap_HR_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/common_SNP/overlap_HR_published.txt"
  titre_intensity_under_HR_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/common_SNP/intensity_under_HR_published.txt"
   titre_f5
  


 titre_graphe_chromosomes <- paste0(dossier_graphes,"rec_chromosomes.png")
 titre_graphe_domestication <-   paste0(dossier_graphes,"domestication.png")
 titre_graphe_HR_partages <-  paste0(dossier_graphes,"HR_shared.png")
 

 titre_domestication <- paste0(dossier_amont,"domestication.csv")
 titre_csre_4Mb <-"/work/adanguy/these/pipeline/020820/csre/map_4mb_csre.txt"
 titre_resume <- "/work/adanguy/these/pipeline/020820/PHASE/PHASE_summary_outputs.txt"
 titre_graphe_domestication <-  "/work/adanguy/these/pipeline/020820/graphs/domestication.png"



 titre_landraces2="/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/supplementary_tab2.txt"
 titre_HAPFLK_tree_SNP <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/hapflk_tree_genome.txt"
 titre_pairwise_FST_matrix_haplotypic_blocks_graph <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/pairwise_FST_matrix_haplotypic_blocks.txt"
 titre_matrice_FST_png <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/pairwise_FST_matrix_haplotypic_blocks.png"
 titre_landraces <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/landraces.txt"
 titre_figure1 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure1.tiff"

 titre_landraces_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/supplementary_tab2.txt"
 titre_HAPFLK_tree_SNP <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/HPAFLK_tree_genome_SNP.txt"
 titre_matrice_FST_blocs_png <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphs/pairwise_FST_matrix_haplotypic_blocks.png"
 titre_figure1 <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/figure1.tiff"
 


head(fread(titre_correlations_4mb_published))
head(fread(titre_FST_published))
head(fread(titre_correlations_mixed_models_published))
head(fread(titre_HR_published))
read.table(titre_matrice_FST_blocs)
head(fread(titre_genetic_maps_published))
head(fread(titre_maps_4Mb))
head(fread(titre_overlap_HR_published))
head(fread(titre_intensity_under_HR_published))



convert_to_adjacency_matrix <- function(dissimilarities, seuil){
  
  # Chnage format
  dissimilarities <- as.matrix(dissimilarities)
  # Trick : dissimilarity > threshold are replaced by 100
  dissimilarities[dissimilarities < seuil] <- -100
  #  dissimality <= threshold are replace by 1
  dissimilarities[dissimilarities >= seuil] <- 1
  # Trick : dissimilarity > threshold are replaced by 0
  dissimilarities[dissimilarities == -100] <- 0
  
  # Sortie : an ajdency matrix (0 and 1)
  return(dissimilarities)
  
}


if (repertoire=="common_SNP"){
  
  
  method="common_SNPs"
} else {
  
  method="population_specific_SNPs"
}


R3="#F8766D"
R2a="#00Bf6D"
C="#00B0F6"
R2b="#A3A500"
R1="#E76BF3"
regions <- c(R1,R2a,C,R2b,R3)


WEEE="#FDE725FF"
EEWA="#B4DE2CFF"
WAEA="#6DCD59FF"
WEWA="#35B779FF"
EEEA="#1F9E89FF"
WEEA="#26828EFF"
populations <- c(WEEE,EEWA,WAEA,WEWA,EEEA,WEEA)


WE="#D7191C"
EE="#FDAE61"
WA="#ABDDA4"
EA="#2B83BA"
populations2 <- c(WE,EE,WA,EA)


inf0 ="#F8766D"
on0="#619CFF"
over0="#00BA38"
slopes_signif=c(inf0, on0, over0)


color_HR=viridis(5)

##### Figure 2 : recombination profiles


 f2a <- fread(titre_maps_4Mb) %>%
   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
   filter(population=="CsRe") %>%
   filter(chr=="3B") %>%
   group_by(population, region) %>%
   slice(n()) %>%
   ungroup() %>%
   mutate(region=lead(region)) %>%
   filter(region !="R1") %>%
   na.omit() %>%
   rbind(., fread(titre_maps_4Mb) %>%
           filter(population=="CsRe") %>%
           filter(chr=="3B")) %>% 
   ggplot(aes(x=posl, y=rec_rate, col=region)) +
   geom_line() +
   ylab("CsRe Bayesian rec rates (cM/Mb)") +
   xlab("physical position (Mb)") +
   ggtitle("Meiotic") +
   theme_light() +
   theme(plot.title = element_text(hjust = 0.5, size=16),
         axis.text.y = element_text(size=12),
         axis.title.y = element_text(size=14),
         axis.text.x = element_text(size=12),
         axis.title.x = element_text(size=14),
         legend.text=element_text(size=12),
         legend.title=element_text(size=14),
         legend.position = "none")+
   scale_colour_manual(values=regions)+
   scale_x_continuous(labels=function(x)x/1e6) 


 # f2b <- fread(titre_maps_4Mb) %>%
 #   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
 #   filter(population!="CsRe") %>%
 #   filter(chr=="3B") %>%
 #   group_by(population, region) %>%
 #   slice(n()) %>%
 #   ungroup() %>%
 #   mutate(region=lead(region)) %>%
 #   filter(region !="R1") %>%
 #   na.omit() %>%
 #   rbind(., fread(titre_maps_4Mb) %>%
 #           filter(population!="CsRe") %>%
 #           filter(chr=="3B")) %>%
 #   ggplot(aes(x=posl, y=rec_rate, col=region)) + 
 #   geom_line() +
 #   ylab(expression(paste(rho," (/kb)"))) +
 #   xlab("physical position (Mb)") +
 #   ggtitle("LD-based") +
 #   theme_light() +
 #   theme(plot.title = element_text(hjust = 0.5)) +
 #   facet_grid(population~.) +
 #   scale_y_continuous(trans="log10")+
 #   theme(plot.title = element_text(hjust = 0.5, size=16),
 #         strip.text.y = element_text(size = 12, color="black"),
 #         axis.text.y = element_text(size=12),
 #         axis.title.y = element_text(size=14),
 #         axis.text.x = element_text(size=12),
 #         axis.title.x = element_text(size=14),
 #         legend.text=element_text(size=12),
 #         legend.title=element_text(size=14),
 #         strip.background = element_rect(color="black", fill="white"))+ 
 #   scale_colour_manual(values=regions)+
 #   scale_x_continuous(labels=function(x)x/1e6)


 # f2b <- fread(titre_maps_4Mb) %>%
 #   mutate(rec_rate=rec_rate*1e3) %>%
 #   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
 #   filter(population!="CsRe") %>%
 #   filter(chr=="3B") %>%
 #   group_by(population, region) %>%
 #   slice(n()) %>%
 #   ungroup() %>%
 #   mutate(region=lead(region)) %>%
 #   filter(region !="R1") %>%
 #   na.omit() %>%
 #   rbind(., fread(titre_maps_4Mb) %>%
 #           mutate(rec_rate=rec_rate*1e3) %>%
 #           filter(population!="CsRe") %>%
 #           filter(chr=="3B")) %>%
 #   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
 #   ggplot(aes(x=posl, y=rec_rate, col=region)) + 
 #   geom_line() +
 #   ylab(expression(paste(rho," (/kb)"))) +
 #   xlab("physical position (Mb)") +
 #   ggtitle("LD-based") +
 #   theme_light() +
 #   theme(plot.title = element_text(hjust = 0.5)) +
 #   facet_grid(population~., scales="free") +
 #   theme(plot.title = element_text(hjust = 0.5, size=16),
 #         strip.text.y = element_text(size = 12, color="black"),
 #         axis.text.y = element_text(size=12),
 #         axis.title.y = element_text(size=14),
 #         axis.text.x = element_text(size=12),
 #         axis.title.x = element_text(size=14),
 #         legend.text=element_text(size=12),
 #         legend.title=element_text(size=14),
 #         panel.grid.minor = element_blank(),
 #         strip.background = element_rect(color="black", fill="white"))+ 
 #   scale_colour_manual(values=regions)+
 #   scale_x_continuous(labels=function(x)x/1e6)+
 #   scale_y_continuous(breaks = c(0,0.01,0.02,0.03,0.04,0.05,0.06), labels = c("0.0","0.1","0.2","0.3","0.4","0.5","0.6"))
 
 
 
 f2b <- fread(titre_maps_4Mb) %>%
   mutate(rec_rate=rec_rate*1e3) %>%
   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
   filter(population!="CsRe") %>%
   filter(chr=="3B") %>%
   group_by(population, region) %>%
   slice(n()) %>%
   ungroup() %>%
   mutate(region=lead(region)) %>%
   filter(region !="R1") %>%
   na.omit() %>%
   rbind(., fread(titre_maps_4Mb) %>%
           mutate(rec_rate=rec_rate*1e3) %>%
           filter(population!="CsRe") %>%
           filter(chr=="3B")) %>%
   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
   ggplot(aes(x=posl, y=rec_rate, col=region)) + 
   geom_line() +
   ylab(expression(paste(rho, " (/kb)"))) +
   xlab("physical position (Mb)") +
   ggtitle("LD-based") +
   theme_light() +
   theme(plot.title = element_text(hjust = 0.5)) +
   facet_grid(population~.) +
   theme(plot.title = element_text(hjust = 0.5, size=16),
         strip.text.y = element_text(size = 12, color="black"),
         axis.text.y = element_text(size=12),
         axis.title.y = element_text(size=14),
         axis.text.x = element_text(size=12),
         axis.title.x = element_text(size=14),
         legend.text=element_text(size=12),
         legend.title=element_text(size=14),
         strip.background = element_rect(color="black", fill="white"))+ 
   scale_colour_manual(values=regions)+
   scale_x_continuous(labels=function(x)x/1e6)+
   scale_y_continuous(trans="log10", breaks = c(1e-4, 1e-3,1e-2), labels = c(expression(10^-4),expression(10^-3),expression(10^-2)) ) 
   

    


f2 <- ggdraw() +
  draw_plot(f2a, x=0, y=0.22, width=0.33, height = 0.565) +
  draw_plot(f2b, x=0.33, y=0,width=0.66, height = 1) 

f2

tiff(titre_f2, units="in", width = 9, height = 6, res=200, compression = "lzw")
f2
dev.off()



##### Figure 3 : scatter plots and boxplots correlations
minimum_WE <- fread(titre_maps_4Mb) %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  pivot_wider(id_cols = c("chr","region","posl","posr"), names_from = "population", values_from = "rec_rate") %>%
  dplyr::select(WE) %>%
  unlist() %>%
  as.vector() %>%
  min(na.rm=T) * 1e3

f3all <- fread(titre_maps_4Mb) %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  pivot_wider(id_cols = c("chr","region","posl","posr"), names_from = "population", values_from = "rec_rate") %>%
  mutate(WE=1e3*WE) %>%
  ggplot(aes(x=CsRe, y=WE, col=region)) +
  geom_point(size=0.5)+
  theme_light() +
  scale_x_continuous(trans="log10", breaks = c(1e-2, 1e-1,1, 10), labels = function(x) ifelse(x == 1, "1", x)) +
  scale_y_continuous(trans="log10", breaks = c(1e-5,1e-3,1e-1), limits = c(minimum_WE, 1e-1), labels =c(expression(10^-5),expression(10^-3),expression(10^-1)))  +
  theme(plot.title = element_text(hjust = 0.5, size=14, vjust=-1),
        axis.text.y =element_text(size=14),
        axis.title.y = element_blank(),
        axis.text.x = element_blank() ,
        axis.title.x = element_blank(),
        legend.text=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.minor = element_blank())+
  coord_fixed(ratio=0.65)+
  ggtitle(("all regions"))


f3all


# R1

f3_cor_regions <- function(nom_region, col_region ) {
  graphe <- fread(titre_maps_4Mb) %>%
    mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
    pivot_wider(id_cols = c("chr","region","posl","posr"), names_from = "population", values_from = "rec_rate") %>%
    mutate(WE=1e3*WE) %>%
    ggplot(aes(x=CsRe, y=WE, col=region)) +
    scale_colour_manual(values=regions)+
    geom_point(data=fread(titre_maps_4Mb) %>%
                 mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
                 pivot_wider(id_cols = c("chr","region","posl","posr"), names_from = "population", values_from = "rec_rate")%>%
                 mutate(WE=1e3*WE) , aes(x=CsRe, y=WE), col="white", alpha=1, size=0.5) +
    geom_point(data=fread(titre_maps_4Mb) %>%
                 mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
                 pivot_wider(id_cols = c("chr","region","posl","posr"), names_from = "population", values_from = "rec_rate") %>%
                 mutate(WE=1e3*WE), aes(x=CsRe, y=WE), col="lightgrey", alpha=1, size=0.5)+
    geom_point(data=fread(titre_maps_4Mb) %>%
                 mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
                 pivot_wider(id_cols = c("chr","region","posl","posr"), names_from = "population", values_from = "rec_rate") %>%
                 mutate(WE=1e3*WE) %>%
                 filter(region==!!nom_region), aes(x=CsRe, y=WE), col=col_region, size=0.5)  +
    scale_x_continuous(trans="log10", breaks = c(1e-2, 1e-1,1, 10), labels = c(expression(10^-2),expression(10^-1),1,10)) +
    
    scale_y_continuous(trans="log10", breaks = c(1e-5,1e-3,1e-1), limits = c(minimum_WE, 1e-1), labels =c(expression(10^-5),expression(10^-3),expression(10^-1)))  +
    theme_light() +
    xlab("") +
    ylab("")  + 
    theme(plot.title = element_text(hjust = 0.5, size=14, vjust=-1),
          axis.text.y =element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank() ,
          axis.title.x = element_blank(),
          legend.text=element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm"),
          panel.grid.minor = element_blank())+  coord_fixed(ratio=0.65)+
    ggtitle((nom_region))
  
  return(graphe)
  
  
}

#scale_x_continuous(trans="log10", breaks = c(1e-2, 1e-1,1, 10), labels = function(x) ifelse(x == 1, "1", x)) +
  
f3R1 <- f3_cor_regions(nom_region = "R1", col_region = R1) 
f3R2a <- f3_cor_regions(nom_region = "R2a", col_region = R2a) + theme(axis.text.y =element_text(size=14))
f3R2b <- f3_cor_regions(nom_region = "R2b", col_region = R2b) + theme(axis.text.y =element_text(size=14)) + theme(axis.text.x =element_text(size=14))
f3R3 <- f3_cor_regions(nom_region = "R3", col_region = R3)   + theme(axis.text.x =element_text(size=14))
f3C <- f3_cor_regions(nom_region = "C", col_region = C)






WE7DR3 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="7D" & region=="R3" & population=="WE") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector() %>%
  round(2)
EE7DR3 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="7D" & region=="R3" & population=="EE") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector()%>%
  round(2)
WA7DR3 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="7D" & region=="R3" & population=="WA") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector()%>%
  round(2)
EA7DR3 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="7D" & region=="R3" & population=="EA") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector()%>%
  round(2)

WE2AR1 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="2A" & region=="R1" & population=="WE") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector() %>%
  round(2)
EE2AR1 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="2A" & region=="R1" & population=="EE") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector()%>%
  round(2)
WA2AR1 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="2A" & region=="R1" & population=="WA") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector()%>%
  round(2)
EA2AR1 <- fread(titre_correlations_4mb_published) %>%
  filter(chr=="2A" & region=="R1" & population=="EA") %>%
  dplyr::select(correlation) %>%
  unlist() %>%
  as.vector()%>%
  round(2)



f3_cor_csre_rho <- fread(titre_correlations_4mb_published) %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
  
  filter(nb_windows_4Mb>=5) %>%
  ggplot(aes(x=population, y=correlation, col=region, shape="test"))+
  geom_boxplot(show.legend = F) +
  theme_light() +
  geom_point(position=position_jitterdodge(jitter.width=0.1)) +
  ylab(expression(paste("correlation ", rho, " and CsRe")))  +
  xlab("") +
  theme(axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.title.x = element_text(size=16),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        legend.key = element_rect(color = NA, fill = NA))+
  scale_colour_manual(values=regions) +
  geom_point(aes(x=1+0.3,y=WE7DR3, shape="7DR3", col="black"), size=3, col="black")+
  geom_point(aes(x=2+0.3,y=EE7DR3, shape="7DR3", col="black"), size=3, col="black")+
  geom_point(aes(x=3+0.3,y=WA7DR3, shape="7DR3", col="black"), size=3, col="black")+
  geom_point(aes(x=4+0.3,y=EA7DR3, shape="7DR3", col="black"), size=3, col="black")+
  geom_point(aes(x=1-0.3,y=WE2AR1, shape="2AR1", col="black"), size=3, col="black")+
  geom_point(aes(x=2-0.3,y=EE2AR1, shape="2AR1", col="black"), size=3, col="black")+
  geom_point(aes(x=3-0.3,y=WA2AR1, shape="2AR1", col="black"), size=3, col="black")+
  geom_point(aes(x=4-0.3,y=EA2AR1, shape="2AR1", col="black"), size=3, col="black")+
  scale_shape_manual(name = "introgressions",
                     breaks = c("none", "2AR1","7DR3"),
                     values = c(8,4, 20) )+
  guides( shape = guide_legend(order = 2),col = guide_legend(order = 1))

f3_cor_csre_rho


f3 <-ggdraw() +
  draw_plot(f3R2b, x=0+0.02,y=0.12, width=0.2, height = 0.2, scale=1.7) +
  draw_plot(f3R2a, x=0+0.02,y=0.45, width=0.2, height = 0.2, scale=1.5) +
  draw_plot(f3all, x=0+0.02,y=0.75, width=0.2, height = 0.2, scale=1.5) +
  draw_plot(f3R3, x=0.15+0.06,y=0.12, width=0.2, height = 0.2, scale=1.7) +
  draw_plot(f3C, x=0.15+0.06,y=0.45, width=0.2, height = 0.2, scale=1.5) +
  draw_plot(f3R1, x=0.15+0.06,y=0.75, width=0.2, height = 0.2, scale=1.5) +
  draw_plot(f3_cor_csre_rho, x = 0.45, width=0.56) +
  draw_label(label=expression(paste("WE ", rho, " (/kb)")), c(0.01), c(0.5), size =14, fontfamily="sans", fontface = "plain", angle=90)+
  draw_label(label=c("CsRe Bayesian recombination rates (cM/Mb)"), c(0.23), c(0.02), size = 14, fontfamily="sans", fontface = "plain")

f3


tiff(titre_f3, units="in", width =12, height =6, res=200, compression = "lzw")
f3
dev.off()

# Figure 4 : correlation per paire of pop






matrice_fst_blocs <- read.table(titre_matrice_FST_blocs) %>% as.matrix() %>%
  melt(value.name="FST") %>%
  rename(P1=Var1,P2=Var2) %>%
  mutate(P=paste0(P1,"-",P2)) %>%
  mutate(P = factor(P, levels=c("WE-EE","EE-WA","WA-EA","WE-WA","EE-EA","WE-EA"))) %>%
  na.omit() %>%
  mutate(FST=100*round(FST,3)) %>%
  arrange(FST) %>%
  dplyr::select(P, FST)

# 

fread(titre_correlations_mixed_models_published) %>% 
  filter(region !="C") %>%
  mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
  mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
  group_by(P1,P2) %>%
  summarise(c=round(mean(cor_lambda), digits=2)) %>%
  arrange(desc(c))


t <- pairwise.t.test(fread(titre_correlations_mixed_models_published) %>% 
                       filter(region !="C") %>%
                       dplyr::select(cor_lambda) %>%
                       unlist() %>%
                       as.vector(),
                     fread(titre_correlations_mixed_models_published) %>%
                       filter(region !="C") %>%
                       mutate(P=paste0(P1,"-",P2)) %>%
                       dplyr::select(P) %>%
                       unlist() %>%
                       as.vector(),
                     p.adj = "bonf")$p.value
m=rbind.data.frame(NA,t)
m=cbind.data.frame(m,NA)
colnames(m) <- c(colnames(m)[1], row.names(m)[-1])
row.names(m) <- c(colnames(m)[1], row.names(m)[-1])
diag(m) <- 1
adjency_matrix<- convert_to_adjacency_matrix(m,0.05)
graphe <- graph_from_adjacency_matrix(adjency_matrix, mode = "lower", weighted = T, diag = F)



tiff(titre_graph_significativity_boxplot_pairs_of_populations, units="in", width = 6, height = 6, res=200, compression = "lzw")
plot(graphe,vertex.size=3, main="Connected Components")
dev.off()
adjency_matrix


if (repertoire =="all"){
  
  texte <- cbind(data.frame(P1=combn(c("WE","EE","WA","EA"), m=2)),data.frame(P1=combn(rev(c("WE","EE","WA","EA")), m=2))) %>%
    t() %>%
    as.data.frame() %>%
    rename(P1=V1,P2=V2) %>%
    filter(P1 !=P2) %>%
    mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
    mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
    arrange(P1,P2)  %>%
    filter((P1=="WE" & P2=="EE") | 
             (P1=="WE" & P2=="WA") |
             (P1=="WE" & P2=="EA") |
             (P1=="EE" & P2=="WA") |
             (P1=="EE" & P2=="EA") |
             (P1=="WA" & P2=="EA")) %>%
    mutate(signif=c("a","b","cd","b","d","c")) %>%
    rowwise() %>%
    mutate(P=paste0(P1,"-",P2)) %>%
    mutate(P = factor(P, levels=c("WE-EE","EE-WA", "WE-WA" ,"WA-EA","EE-EA", "WE-EA" ))) %>%
    arrange(P)%>%
    dplyr::select(P, signif)
  
  
} else if (repertoire=="common_SNP"){
  
  
  texte <- cbind(data.frame(P1=combn(c("WE","EE","WA","EA"), m=2)),data.frame(P1=combn(rev(c("WE","EE","WA","EA")), m=2))) %>%
    t() %>%
    as.data.frame() %>%
    rename(P1=V1,P2=V2) %>%
    filter(P1 !=P2) %>%
    mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
    mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
    arrange(P1,P2)  %>%
    filter((P1=="WE" & P2=="EE") | 
             (P1=="WE" & P2=="WA") |
             (P1=="WE" & P2=="EA") |
             (P1=="EE" & P2=="WA") |
             (P1=="EE" & P2=="EA") |
             (P1=="WA" & P2=="EA")) %>%
    mutate(signif=c("a","ab","c","b","c","c")) %>%
    rowwise() %>%
    mutate(P=paste0(P1,"-",P2)) %>%
    mutate(P = factor(P, levels=c("WE-EE","EE-WA", "WE-WA" ,"WA-EA","EE-EA", "WE-EA" ))) %>%
    arrange(signif)%>%
    dplyr::select(P, signif)
  
  
}




 
 f4 <- fread(titre_correlations_mixed_models_published) %>% 
   filter(method==method) %>%
   mutate(P=paste0(P1,"-",P2)) %>%
   mutate(P = factor(P, levels=c("WE-EE","WE-WA", "EE-WA","WA-EA","EE-EA","WE-EA"))) %>%
   ggplot(aes(x=P,y=cor_lambda, fill=P)) + 
   geom_boxplot()+
   scale_fill_manual(values=c(WEEE,EEWA,WAEA,WEWA,EEEA,WEEA), 
                     labels=matrice_fst_blocs$FST,name="")+
   theme_light()+
   theme(plot.title = element_text(hjust = 0.5, size=24),
         axis.text.y = element_text(size=22),
         axis.title.y = element_text(size=24),
         axis.text.x = element_text(size=22, angle = 45, hjust=1),
         axis.title.x = element_text(size=24),
         legend.text=element_text(size=22, colour="red"),
         legend.title=element_text(size=24, colour="red"),
         legend.key = element_rect(fill = "white"),
         legend.key.width = unit(3, "cm"),
         legend.key.height = unit(2, "cm"),
         legend.position = "bottom",
         legend.box = "horizontal",
         legend.direction = "horizontal", 
         legend.spacing.x = unit(1.7, 'cm'))+
   guides(fill = guide_legend(label.hjust = 0.3, 
                              label.position = "top",
                              title.position = "top",
                              title.hjust = 0.5,
                              nrow = 1,
                              override.aes=list(fill=NA, col=NA, shape=NA))) +
   
   xlab("") +
   stat_summary(geom = 'text', 
                label = texte$signif, 
                fun.y = max, 
                vjust = -1, 
                size=10)  +
   scale_y_continuous(limits = c(min(fread(titre_correlations_mixed_models_published)$cor_lambda),
                                 max(fread(titre_correlations_mixed_models_published)$cor_lambda)+0.1))+
   ylab(expression(paste("correlation log10(", lambda, ")")))



 arrow_p <- 
   ggplot(tips) +
   geom_segment(aes(x = 0, y = 0, xend = 42, yend = 0),
                col = "red",
                arrow = arrow(length = unit(0.3, "cm"))) +
   coord_cartesian(ylim = c(0, 0.75), xlim = c(0,50)) +
   theme_void()
 
 
 f4f <-ggdraw() +
   draw_plot(f4) +
   draw_plot(arrow_p, x = 0.1, y = 0.05)+
   draw_plot_label(c("Fst (%)"), c(-0.04), c(0.13), size = 25, color="red", fontface = "plain", family = "sans")
 f4f

# f4f <- fread(titre_correlations_mixed_models_published) %>% 
#   filter(method==method) %>%
#   filter(region !="C") %>%
#   mutate(P=paste0(P1,"-",P2)) %>%
#   full_join(matrice_fst_blocs, by="P") %>%
#   mutate(P = factor(P, levels=c("WE-EE","EE-WA", "WE-WA" ,"WA-EA","EE-EA", "WE-EA" ))) %>%
#   mutate(FST =as.factor(FST)) %>%
#   ggplot(aes(x=FST,y=cor_lambda, fill=P, group=P)) + 
#   geom_boxplot()+
#   scale_fill_manual(values=c(WEEE,EEWA,WAEA,WEWA,EEEA,WEEA),name="")+
#   theme_light()+
#   theme(plot.title = element_text(hjust = 0.5, size=24),
#         axis.text.y = element_text(size=22),
#         axis.title.y = element_text(size=24),
#         axis.text.x = element_text(size=22),
#         axis.title.x = element_text(size=24),
#       legend.key.size = unit(1, "cm"),
#       legend.position = "bottom",
#       legend.spacing.x = unit(0.5, 'cm'),
#         legend.text=element_text(size=22))+
#   
#   xlab("") +
#   guides(fill = guide_legend(nrow=1))+
#   stat_summary(geom = 'text', 
#                label = texte$signif, 
#                fun.y = max, 
#                vjust = -1, 
#                size=10)  +
#   scale_y_continuous(limits = c(min(fread(titre_correlations_mixed_models_published) %>% filter(region != "C") %>%   filter(method=="population_specific_SNPs") %>%dplyr::select(cor_lambda)),
#                                 max(fread(titre_correlations_mixed_models_published) %>% filter(region != "C")  %>% filter(method=="population_specific_SNPs") %>% dplyr::select(cor_lambda))+0.1))+
#   ylab(expression(paste("correlation log10(", lambda, ")")))
# f4f



f4f
tiff(titre_f4, units="in", width = 12, height = 9, res=200, compression = "lzw")
f4f
dev.off()


# Figure 5 : HR


f5a <- fread(titre_overlap_HR_published) %>%
  mutate(groupe=case_when(populations%in%c("WE","EE","WA","EA")~"population specific",
                          populations %in% c( "WEEE","WEWA","WEEA","EEWA","EEEA","WAEA") ~ "2 populations +" ,
                          populations %in% c("WEEEWA","WEEEEA","WEWAEA","EEWAEA")~"3 populations +",
                          populations %in% c("WEEEWAEA")~ "4 populations"))%>%
  mutate(groupe=factor(groupe, levels = rev(c("4 populations","3 populations +","2 populations +","population specific")))) %>%
  group_by(iteration) %>%
  mutate(nbHRIs_ref=ifelse(iteration=="ref", nbHRIs, NA)) %>%
  mutate(nb_total_cliques_ref=ifelse(iteration=="ref", nb_total_cliques, NA)) %>%
  ungroup() %>%
  group_by(populations) %>%
  mutate(nbHRIs_ref=max(nbHRIs_ref, na.rm=T)) %>%
  mutate(nb_total_cliques_ref=max(nb_total_cliques_ref, na.rm=T)) %>%
  filter(iteration != "ref") %>%
  group_by(groupe, iteration, nb_total_cliques_ref,nb_total_cliques) %>%
  summarise(vref=sum(nbHRIs_ref), v=sum(nbHRIs)) %>%
  arrange(iteration, desc(groupe)) %>%
  mutate(groupe2=ifelse(groupe %in% "population specific", "population specific", "2ormore"))%>%
  group_by(iteration, groupe2) %>%
  mutate(vref=cumsum(vref), v=cumsum(v)) %>%
  ungroup() %>%
  mutate(vref=vref/nb_total_cliques_ref, v=v/nb_total_cliques) %>%
  dplyr::select(groupe, iteration, vref, v) %>%
  ggplot(aes(x=groupe, y =v, fill="temp")) + geom_boxplot(col="grey", show.legend = T)  +
  scale_fill_manual(values="grey", name="", labels="no sharing")+
  geom_point(aes(x=groupe, y=vref, col=groupe), shape=18, size=5, show.legend = F) +
  scale_color_manual(values=c("4 populations"=color_HR[4],
                              "3 populations +"=color_HR[3],
                              "2 populations +"=color_HR[2],
                              "population specific"=color_HR[1]                              ), guide=F) +
  ylab("Percentage of HRIs shared by") +
  facet_grid(.~groupe, scales="free")+
  xlab("")  +
  theme_light() +
  geom_text(aes(x=groupe, y=vref, col=groupe, label=round(vref, digits=2)),hjust=0.8, vjust=-1, check_overlap = TRUE, size=10) +
  theme(legend.position = "right",
        strip.text.x = element_text(size = 26, color="black"),
        strip.background = element_rect(color="black", fill="white"),
        axis.title.y = element_text(size=26),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=24),
        axis.title.x = element_blank(),
        legend.text=element_text(size=26),
        legend.title=element_text(size=22),
        panel.grid.minor = element_blank())
f5a






f5b <- fread(titre_intensity_under_HR_published) %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
  mutate(population_ref=factor(population_ref, levels=c("WE","EE","WA","EA"))) %>%
  group_by(population, population_ref, pos) %>%
  summarise(rel_lambda=median(rel_lambda)) %>%
  ggplot(aes(x=pos, y=rel_lambda, col=population)) +
  geom_line(size=2) +
  facet_wrap(population_ref~., ncol=4, scales="free") +
  theme_light() +
  xlab("distance to HRI center (kb)") +
  ylab(expression(paste(lambda, " tested pop / ",lambda," reference pop (%)"))) +
  scale_colour_manual(values=populations2, name="tested\npopulation") +
  scale_x_continuous(breaks = seq(-1e5, 1e5, 1e5), labels=function(x)x/1e3)  +
  scale_y_continuous(labels=function(x)x*100, limits=c(0.11,1), trans="sqrt", breaks = c(0,0.25,0.5,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=24),
        plot.margin = unit(c(0,0,1,0), "cm"),
        strip.text.x = element_text(size = 26, color="black"),
        strip.background = element_rect(color="black", fill="white"),
        axis.title.y = element_text(size=26),
        axis.title.x = element_text(size=26),
        axis.text.y = element_text(size=24),
        legend.text=element_text(size=26),
        legend.title=element_text(size=26),
        panel.grid.minor = element_blank()) 



f5b


f5 <- ggarrange(f5a,f5b, nrow=2, labels = c("A","B"), font.label = list(size=30))
f5


tiff(titre_f5, units="in", width = 18, height = 18, res=200, compression = "lzw")
f5
dev.off()



# Figure 6 : relation correlation lambda dn FST per genomic regin 1AR1...7DR3




mod_lambda_haplotypic_blocks = lmList( cor_lambda ~ FST | chrregion, data=fread(titre_correlations_mixed_models_published) %>%
                                         inner_join(fread(titre_FST_published), by=c("chr"="chr","region"="region","P1"="P1","P2"="P2")) %>%
                                         filter(method=="haplotypic_blocks") %>% mutate(chrregion=paste0(chr,region)) %>%
                                         na.omit() %>%
                                         filter(region !="C"))


slopes_lambda_haplotypic_blocks=data.frame(intervals(mod_lambda_haplotypic_blocks)[,,2]) %>%
  rename(slope=est.) %>%
  rownames_to_column(var="chrregion") %>%
  mutate(chr=substr(chrregion,1,2)) %>%
  mutate(region=substr(chrregion,3,5)) %>%
  mutate( confidence_interval = case_when(upper<0~"<0",
                                          lower>0~">0",
                                          upper>0 & lower<0 ~ "include 0")) %>%
  mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0"))) 







f6a <- fread(titre_correlations_mixed_models_published) %>%
  inner_join(fread(titre_FST_published), by=c("chr"="chr","region"="region", "P1"="P1","P2"="P2")) %>%
  filter(method=="haplotypic_blocks") %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  filter(region!="C") %>%
  inner_join(slopes_lambda_haplotypic_blocks, by=c("chr"="chr","region"="region"))%>%
  ggplot(aes(x=FST, y=cor_lambda)) + 
  geom_point() +
  geom_smooth(method="lm", se=F, aes(col=confidence_interval)) +
  scale_color_manual(values=slopes_signif) +
  facet_grid(region~chr, scales="free") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10))+
  theme( axis.text.y = element_text(size=10)) +
  theme(legend.position="bottom") +
  ylab(expression(paste("correlation of log10(", lambda,")"))) +
  xlab("\nFst (%, from haplotypic alleles)") +
  guides(colour=guide_legend(title="95% confidence interval of slope estimate : ", override.aes = list(size = 3)))+
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+
  scale_x_continuous( labels=function(x)x*100,breaks=seq(0,1,0.05))+
  theme(strip.text.y = element_text(size = 22, color="black"),
        strip.text.x = element_text(size = 22, color="black"),
        strip.background = element_rect(color="black", fill="white"),
        axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        legend.text=element_text(size=22),
        legend.title=element_text(size=22),
        panel.grid.minor = element_blank())
f6a





f6b <- slopes_lambda_haplotypic_blocks %>%
  mutate( cr = fct_reorder(chrregion, slope)) %>%
  ggplot(aes(x=slope,y=cr,colour=confidence_interval)) +
  scale_color_manual(values=slopes_signif) +
  geom_segment(aes(x=lower,xend=upper,yend=chrregion),colour='gray', size=1) +
  geom_point() +
  geom_vline(xintercept=0,col=slopes_signif[2], size=1) +
  xlab(label='Fst effect') + ylab( label = "genomic region") +
  theme_light()+
  theme( axis.text.y = element_text(size=10)) +
  theme( axis.text.x = element_text(size=18)) +
  theme(legend.position="none") +
  guides(colour=guide_legend("95% confidence interval of slope estimate: "))+
  theme(axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        legend.text=element_text(size=22),
        legend.title=element_text(size=22),
        panel.grid.minor = element_blank())




f6 <- ggdraw() +
  draw_plot(f6a, x=-0.33, y=0,1, 0.5, width =1.35) +
  draw_plot(f6b,x=0.7,y=0,1,width=0.3)+
  draw_plot_label(c("A", "B"), c(0,0.68), c(0.85,1), size = 30)


f6


tiff(titre_f6, units="in", width = 18, height = 9, res=200, compression = "lzw")
f6
dev.off()



# comparison with common SNP


if (repertoire !="common_SNP"){
  
  mod_lambda_haplotypic_blocks_common_SNP = lmList( cor_lambda ~ FST | chrregion, data=fread(titre_correlations_mixed_models_published_common_SNP) %>%
                                                      inner_join(fread(titre_FST_published), by=c("chr"="chr","region"="region","P1"="P1","P2"="P2")) %>%
                                                      filter(method=="haplotypic_blocks") %>% mutate(chrregion=paste0(chr,region)) %>%
                                                      na.omit() %>%
                                                      filter(region !="C"))
  
  slopes_SNP_specific_or_commona <-data.frame(intervals(mod_lambda_haplotypic_blocks_common_SNP)[,,2]) %>%
    rename(slope=est.) %>%
    rownames_to_column(var="chrregion") %>%
    mutate(chr=substr(chrregion,1,2)) %>%
    mutate(region=substr(chrregion,3,5)) %>%
    mutate( confidence_interval = case_when(upper<0~"<0",
                                            lower>0~">0",
                                            upper>0 & lower<0 ~ "include 0")) %>%
    mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0"))) %>%
    mutate(method="common") %>%
    inner_join(slopes_lambda_haplotypic_blocks %>%
                 mutate(method="specific"), by=c("chr"="chr", "region"="region", "chrregion"="chrregion"), suffix=c(".common",".specific")) %>%
    mutate(overlapp_of_conf_interval=ifelse(!(upper.specific <= lower.common | lower.specific >= upper.common), T, F)) %>%
    ggplot(aes(x=slope.specific, y=slope.common)) +
    geom_point() +
    geom_smooth(method = "lm", col="red") +
    geom_abline(slope=1, intercept = 0) +
    theme_light() +
    xlab("Fst effect computed using specific SNP dataset") +
    ylab("Fst effect computed using common SNP dataset")+
    geom_hline(yintercept = 0, col=slopes_signif[2], size=2) +
    geom_vline(xintercept=0, col=slopes_signif[2], size=2) +
    theme(legend.position = "bottom",
          axis.title.y = element_text(size=14, colour="#C77CFF"),
          axis.title.x = element_text( size=14,colour="#F1A340"),
          axis.text.y = element_text( size=14,colour="#C77CFF"),
          axis.text.x = element_text( size=14,colour="#F1A340")) +
    scale_x_continuous(limits = c(-10,2.5))+
    scale_y_continuous(limits = c(-10,2.5))
  
  
  
  
  slopes_SNP_specific_or_commonb <- data.frame(intervals(mod_lambda_haplotypic_blocks_common_SNP)[,,2]) %>%
    rename(slope=est.) %>%
    rownames_to_column(var="chrregion") %>%
    mutate(chr=substr(chrregion,1,2)) %>%
    mutate(region=substr(chrregion,3,5)) %>%
    mutate( confidence_interval = case_when(upper<0~"<0",
                                            lower>0~">0",
                                            upper>0 & lower<0 ~ "include 0")) %>%
    mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0"))) %>%
    mutate(method="common") %>%
    rbind(., slopes_lambda_haplotypic_blocks %>%
            mutate(method="specific")) %>%
    mutate(chrregion=factor(chrregion, levels=slopes_lambda_haplotypic_blocks %>% arrange(slope) %>%
                              dplyr::select(chrregion) %>%
                              unlist() %>%
                              as.vector())) %>%
    ggplot(aes(y = slope, ymin = lower,
               x = chrregion, ymax = upper, group = chrregion, col=method)) +
    geom_point(position = position_dodge(1), size=4) +
    geom_linerange(position = position_dodge(1)) +
    geom_abline(slope=0, intercept = 0, col=slopes_signif[2], size=2) +
    coord_flip()   +
    ylab(label='Fst effect') +
    xlab( label = "genomic region") +
    scale_colour_manual(values=c("#C77CFF","#F1A340"))+
    guides(colour=guide_legend("slope estimate: "))+
    theme_light() +
    theme(axis.title.y = element_text(size=14),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=6),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14),
          panel.grid.minor = element_blank(),
          legend.position = "bottom") 
  
  
  
  slopes_SNP_specific_or_common <- ggdraw() +
    draw_plot(slopes_SNP_specific_or_commonb, x=0.02, y=0.5, height = 0.5, width = 0.95) +
    draw_plot(slopes_SNP_specific_or_commona, x=0.25, y=0., height = 0.5, width = 0.5) +
    draw_plot_label(c("A", "B"), c(0,0), c(1,0.5), size = 20)
  
  titre_slopes_SNP_specific_or_common
  tiff(titre_slopes_SNP_specific_or_common, units="in", width = 12, height = 12, res=200, compression = "lzw")
  slopes_SNP_specific_or_common
  dev.off()
  
  
  
  
  # landraces


  
  # So WE ~ NWE = STRmemb4_3
  # EE ~ SEE = STRmemb4_4
  # WA ~ CAA = STRmemb4_1
  # EA ~ SEA = STRmemb4_2
  
  
  

  # STRUCTURE
  
  WEK4 <- fread(titre_landraces_published) %>%
    filter(population=="WE")%>%
    dplyr::select(ID_LINE, NWE, SEE, CAA, SEA) %>%
    mutate(value_dom=SEE+CAA+SEA) %>%
    mutate(value_dom2=NWE) %>%
    gather(memb, value, -ID_LINE, -value_dom, -value_dom2) %>% 
    rename(STRgp4=memb) %>%
    mutate(STRgp4=factor(STRgp4, levels=c("NWE","SEE","CAA","SEA"))) %>%
    arrange(value_dom, desc(value_dom2)) %>%
    mutate(ID_LINE = factor(ID_LINE, levels=unique(ID_LINE))) %>%
    mutate(pos=as.factor(as.numeric(ID_LINE))) %>%
    ggplot(aes(x=pos, y=value, fill=STRgp4, col=STRgp4) ) + 
    geom_bar(stat="identity", position="fill") +
    theme_void ()+
    scale_colour_manual("K=4", values=populations2) +
    scale_fill_manual("K=4",  values=populations2) +
    theme(legend.key.size = unit(0.5,"line"),
          text=element_text(family="serif", face="bold"),
          legend.box.spacing=unit(1, 'cm'),
          legend.title=element_text(size=15),
          legend.text=element_text(size=10),
          plot.margin = unit(c(1,0.05,0,0.05), "cm"))+
    scale_y_continuous(limits=c(0,0.95),oob = rescale_none) 
  
  
  EEK4 <- fread(titre_landraces_published) %>%
    filter(population=="EE")%>%
    dplyr::select(ID_LINE, NWE, SEE, CAA, SEA) %>%
    mutate(value_dom=NWE) %>%
    mutate(value_dom2=CAA) %>%
    gather(memb, value, -ID_LINE, -value_dom, -value_dom2) %>% 
    rename(STRgp4=memb) %>%
    mutate(STRgp4=factor(STRgp4, levels=c("NWE","SEE","CAA","SEA"))) %>%
    arrange(desc(value_dom), value_dom2) %>%
    mutate(ID_LINE = factor(ID_LINE, levels=unique(ID_LINE))) %>%
    mutate(pos=as.factor(as.numeric(ID_LINE))) %>%
    ggplot(aes(x=pos, y=value, fill=STRgp4, col=STRgp4) ) + 
    geom_bar(stat="identity", position="fill") +
    theme_void ()+
    scale_colour_manual("K=4", values=populations2) +
    scale_fill_manual("K=4",  values=populations2) +
    theme(legend.key.size = unit(0.5,"line"),
          text=element_text(family="serif", face="bold"),
          legend.box.spacing=unit(1, 'cm'),
          legend.title=element_text(size=15),
          legend.text=element_text(size=10),
          plot.margin = unit(c(1,0.05,0,0.05), "cm"))+
    scale_y_continuous(limits=c(0,0.95),oob = rescale_none) 
  
  
  WAK4 <- fread(titre_landraces_published)%>% filter(population=="WA")%>%
    dplyr::select(ID_LINE, NWE, SEE, CAA, SEA) %>%
    mutate(value_dom=SEE) %>%
    mutate(value_dom2=SEA) %>%
    gather(memb, value, -ID_LINE, -value_dom, -value_dom2) %>% 
    rename(STRgp4=memb) %>%
    mutate(STRgp4=factor(STRgp4, levels=c("NWE","SEE","CAA","SEA"))) %>%
    arrange(desc(value_dom), value_dom2) %>%
    mutate(ID_LINE = factor(ID_LINE, levels=unique(ID_LINE))) %>%
    mutate(pos=as.factor(as.numeric(ID_LINE))) %>%
    ggplot(aes(x=pos, y=value, fill=STRgp4, col=STRgp4) ) + 
    geom_bar(stat="identity", position="fill") +
    theme_void ()+
    scale_colour_manual("K=4", values=populations2) +
    scale_fill_manual("K=4",  values=populations2) +
    theme(legend.key.size = unit(0.5,"line"),
          text=element_text(family="serif", face="bold"),
          legend.box.spacing=unit(1, 'cm'),
          legend.title=element_text(size=15),
          legend.text=element_text(size=10),
          plot.margin = unit(c(1,0.05,0,0.05), "cm"))+
    scale_y_continuous(limits=c(0,0.95),oob = rescale_none) 
  
  
  EAK4 <- fread(titre_landraces_published)%>% filter(population=="EA")%>%
    dplyr::select(ID_LINE, NWE, SEE, CAA, SEA) %>%
    mutate(value_dom=SEE+NWE+CAA) %>%
    mutate(value_dom2=CAA) %>%
    gather(memb, value, -ID_LINE, -value_dom, -value_dom2) %>% 
    rename(STRgp4=memb) %>%
    mutate(STRgp4=factor(STRgp4, levels=c("NWE","SEE","CAA","SEA"))) %>%
    arrange(desc(value_dom), value_dom2) %>%
    mutate(ID_LINE = factor(ID_LINE, levels=unique(ID_LINE))) %>%
    mutate(pos=as.factor(as.numeric(ID_LINE))) %>%
    ggplot(aes(x=pos, y=value, fill=STRgp4, col=STRgp4) ) + 
    geom_bar(stat="identity", position="fill") +
    theme_void ()+
    scale_colour_manual("K=4", values=populations2) +
    scale_fill_manual("K=4",  values=populations2) +
    theme(legend.key.size = unit(0.5,"line"),
          text=element_text(family="serif", face="bold"),
          legend.box.spacing=unit(1, 'cm'),
          legend.title=element_text(size=15),
          legend.text=element_text(size=10),
          plot.margin = unit(c(1,0.05,0,0.05), "cm"))+
    scale_y_continuous(limits=c(0,0.95),oob = rescale_none) 
  
  
  K4 <- ggarrange(WEK4, EEK4, WAK4, EAK4, ncol=4, common.legend = T, legend = "none")
  
  
  # tree
  
  d = data.frame(node=1:7, color=c(populations2[3],populations2[4],populations2[2],populations2[1],"black","black","black"))
  
  tree <- read.tree(titre_HAPFLK_tree_SNP) %>%
    ggtree( size=2) %<+% d + aes(color=I(color)) +
    theme_tree() +
    theme(legend.position="none") +
    scale_y_reverse() +
    scale_x_reverse() +
    coord_flip() +
    geom_tiplab(geom="label", offset=-0.005, size=10, align=F, hjust = 0.5, fontface="bold") +
    theme(plot.margin = unit(c(0,0,0,0), "cm")) 
  
  tree <- ggtree::rotate(tree, 6) 

  
  # merge
  
  
  arbre <- ggdraw() +
    draw_plot(tree, x=0, .48, 0.9, .5) +
    draw_plot(K4, 0, 0, 1, .5) +
    draw_plot_label(c("Population tree", "Population structure"), c(-0.09, -0.12), c(1, 0.5), size = 20)
  
  
  
  
  ### FST matrix
  
  dev.off()
  png(titre_matrice_FST_blocs_png, units="in", width = 24, height = 24, res=200)
  par(xpd = TRUE)
  corrplot(round(as.matrix(read.table(titre_matrice_FST_blocs))*100,digits=2),
           method="color", col=c(rep("gray94", times=11000),
                                 rep(populations[1], times=2000),
                                 rep(populations[2], times=1100),
                                 rep(populations[3], times=700),
                                 rep(populations[4], times=5200),
                                 rep(populations[5], times=100),
                                 rep(populations[6], times=1)),  
           type="upper", 
           addCoef.col = "black", 
           tl.col=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA"), 
           tl.srt=45, 
           sig.level = 0.01, 
           insig = "blank", 
           diag=T , is.corr = F, shade.lwd=0,number.digits=1, cl.pos = "n",
           tl.cex=10, 
           number.cex=10,
           mar = c(0, 0, 0, 0))
  dev.off()
  p <- readPNG(titre_matrice_FST_blocs_png)
  
  
  graphe <- ggdraw() +
    draw_image(titre_matrice_FST_blocs_png, 0.7, 0.2, 0.3, 1) +
    draw_plot(arbre, x=0, 0, 0.7, 1) +
    draw_plot_label(c("Fst matrix (%)", ""), c(0.7, 0.7), c(1, 1), size = 20)
  
  tiff(titre_figure1, units="in", width = 15, height = 9, res=200, compression = "lzw")
  graphe
  dev.off()
  
  
  
}




sessionInfo()