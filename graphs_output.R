

#### Graphics

Sys.time()
cat("\n\nGraphs_output.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)



suppressPackageStartupMessages(library(nlme))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))


variables <- commandArgs(trailingOnly=TRUE)

cat("\n Graphics and stats\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

repertoire <- variables[1]
titre_correlations_4mb_published <- variables[2]
titre_FST_published <- variables[3]
titre_correlations_mixed_models_published <- variables[4]
titre_HR_published <- variables[5]
titre_cliques_published<- variables[6]
titre_matrice_fst_blocs<- variables[7]
titre_genetic_maps_published<- variables[8]
titre_maps_4Mb <- variables[9]
titre_f2 <- variables[10]
titre_f3 <- variables[11]
titre_f4 <- variables[12]
titre_f6 <- variables[13]
titre_graph_significativity_boxplot_pairs_of_populations <- variables[14]


# titre_correlations_4mb_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/correlations_4mb_published.txt"
# titre_FST_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/FST_published.txt"
# titre_correlations_mixed_models_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/correlations_mixed_models_published.txt"
# titre_HR_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/HR_published.txt"
# titre_cliques_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/titre_cliques_published.txt"
# titre_matrice_fst_blocs <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/pairwise_FST_matrix_haplotypic_blocks.txt"
# titre_genetic_maps_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/genetic_maps_published.txt"
# titre_maps_4Mb <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/historical_and_meiotic_maps_4mb_published.txt"
# titre_f3 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure3.tiff"
# titre_f2 <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure2.tiff"
# titre_f4 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure4.tiff"
# titre_f6 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure6.tiff"
# 
# 



# titre_resume <- paste0(dossier,"PHASE_summary_outputs.txt")
# titre_csre_4Mb <- paste0(dossier,"map_4mb_csre.txt")
# titre_landraces_4Mb <- paste0(dossier,"map_4mb_landraces.txt")
# titre_FST_SNP <-paste0(dossier,"chrregion_FST_SNP.txt")
# titre_FST_blocs <- paste0(dossier,"chrregion_FST_haplotypic_blocks.txt")
# titre_rho_aleatoires <- paste0(dossier,"chrregion_rho_random.txt")
# titre_lambda_aleatoires <-paste0(dossier,"chrregion_lambda_random.txt")
# titre_hrec <- paste0(dossier,"HR.txt")
# titre_csre <- paste0(dossier,"csre_genetic_map.txt")
# titre_domestication <- paste0(dossier,"domestication.csv")
# titre_cliques_HR <- paste0(dossier,"cliques_HR.txt")
# titre_chr_tab <-paste0(dossier,"Decoupage_chr_ble.tab")
# titre_matrice_fst_blocs <- paste0(dossier, "pairwise_FST_matrix_haplotypic_blocks.txt")
# 
# 
# titre_graphe_chromosomes <- paste0(dossier_graphes,"rec_chromosomes.png")
# titre_graphe_domestication <-   paste0(dossier_graphes,"domestication.png")
# titre_graphe_HR_partages <-  paste0(dossier_graphes,"HR_shared.png")
# 

# titre_domestication <- paste0(dossier_amont,"domestication.csv")
# titre_csre_4Mb <-"/work/adanguy/these/pipeline/020820/csre/map_4mb_csre.txt"
# titre_resume <- "/work/adanguy/these/pipeline/020820/PHASE/PHASE_summary_outputs.txt"
# titre_graphe_domestication <-  "/work/adanguy/these/pipeline/020820/graphs/domestication.png"


head(fread(titre_correlations_4mb_published))
head(fread(titre_FST_published))
head(fread(titre_correlations_mixed_models_published))
head(fread(titre_HR_published))
read.table(titre_matrice_fst_blocs)
head(fread(titre_genetic_maps_published))
head(fread(titre_maps_4Mb))



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


inf0 ="#F8766D"
on0="#619CFF"
over0="#00BA38"
slopes_signif=c(inf0, on0, over0)

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


f2b <- fread(titre_maps_4Mb) %>%
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
          filter(population!="CsRe") %>%
          filter(chr=="3B")) %>%
  ggplot(aes(x=posl, y=rec_rate, col=region)) + 
  geom_line() +
  ylab(expression(paste(rho," (/kb)"))) +
  xlab("physical position (Mb)") +
  ggtitle("LD-based") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(population~.) +
  scale_y_continuous(trans="log10")+
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
  scale_x_continuous(labels=function(x)x/1e6)




f2 <- ggarrange(f2a,f2b, ncol=2,nrow=1, widths = c(0.33,0.66), common.legend = T, legend = "bottom")
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
  scale_y_continuous(trans="log10", breaks = c(1e-5,1e-3,1e-1), limits = c(minimum_WE, 1e-1))  +
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
    scale_x_continuous(trans="log10", breaks = c(1e-2, 1e-1,1, 10), labels = function(x) ifelse(x == 1, "1", x)) +
    scale_y_continuous(trans="log10", breaks = c(1e-5,1e-3,1e-1), limits = c(minimum_WE, 1e-1))  +
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
# col <- viridis(7)
# 
# lal <- lal %>% filter(P1 != P2) %>%
#   mutate(chr=substr(chrregion,1,2)) %>%
#   mutate(region=substr(chrregion,3,5)) %>%
#   rename(correlation=estimateur) %>%
#   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
#   filter(region !="C") %>%
#   mutate(p=paste0(P1,"-",P2)) %>% 
#   mutate(p = factor(p, levels=c("WE-EE","EE-WA","WA-EA","WE-WA","EE-EA","WE-EA"))) %>%
#   droplevels() %>%
#   na.omit() %>%
#   dplyr::select(chrregion, chr, region, P1, P2, p, correlation, nb_intervalles, BIC) %>%
#   mutate(color=case_when(P1=="WE" & P2=="EE" ~ col[2],
#                          P1=="WE" & P2 =="WA" ~ col[3],
#                          P1=="WE" & P2 =="EA" ~ col[6],
#                          P1=="EE" & P2 =="WA" ~ col[4],
#                          P1=="EE" & P2 =="EA" ~ col[7],
#                          P1=="WA" & P2 =="EA" ~ col[5])) 
# 
# 
# order <- lal %>% group_by(p) %>%
#   summarise(m=median(correlation)) %>%
#   arrange(desc(m)) %>%
#   dplyr::select(p) %>%
#   unlist() %>% 
#   as.vector()
# 
# lal <- lal %>%
#   mutate(p=factor(p, levels=order))
# 
# lal2 <- lal %>%
#   mutate(P1.1=P2, P2.1=P1) %>%
#   mutate(P1=P1.1, P2=P2.1) %>%
#   dplyr::select(chrregion, chr, region, P1, P2, p, correlation, nb_intervalles, BIC, color) %>%
#   rbind(., lal) %>%
#   mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
#   mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
#   arrange(P1,P2) %>%
#   filter((P1=="WE" & P2=="EE") | 
#            (P1=="WE" & P2=="WA") |
#            (P1=="WE" & P2=="EA") |
#            (P1=="EE" & P2=="WA") |
#            (P1=="EE" & P2=="EA") |
#            (P1=="WA" & P2=="EA"))
# 
# 
# 
# matrice_fst2 <- matrix(unlist(matrice_fst), nrow = 4, ncol = 4, byrow = TRUE,
#                        dimnames = list(colnames(matrice_fst),
#                                        colnames(matrice_fst)))
# 
# matrice_fst3 <- melt(matrice_fst2) %>%
#   rename(P1=Var1,P2=Var2, FST=value) %>%
#   inner_join(lal2, by=c("P1"="P1","P2"="P2")) %>%
#   mutate(p=factor(p, levels=order)) %>%
#   arrange(FST)
# 
# 
# cor_lambda_fst_blocs <-  lal %>%
#   inner_join(FSTb, by=c("chrregion","P1","P2")) %>%
#   na.omit()
# 
# 
# 
# mod = lmList( correlation ~ FST | chrregion, data=cor_lambda_fst_blocs)
# eff.FST.cor_lambda_fst_blocs=data.frame(intervals(mod)[,,2])
# colnames(eff.FST.cor_lambda_fst_blocs)[2]='estimate'
# eff.FST.cor_lambda_fst_blocs$chrregion = rownames(eff.FST.cor_lambda_fst_blocs)
# eff.FST.cor_lambda_fst_blocs <- eff.FST.cor_lambda_fst_blocs %>%
#   mutate( confidence_interval = case_when(upper<0~"<0",
#                                           lower>0~">0",
#                                           upper>0 & lower<0 ~ "include 0")) %>%
#   mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0")))
# 
# 
# cor_lambda_fst_blocs <- cor_lambda_fst_blocs %>% inner_join(eff.FST.cor_lambda_fst_blocs %>% dplyr::select(chrregion, confidence_interval), by="chrregion")
# 
# 
# 
# 
# cor_lambda_fst_snp <- lal %>%
#   inner_join(FSTs, by=c("chrregion","P1","P2")) %>%
#   na.omit()
# 
# 
# mod = lmList( correlation ~ FST | chrregion, data=cor_lambda_fst_snp)
# eff.FST.cor_lambda_snp=data.frame(intervals(mod)[,,2])
# colnames(eff.FST.cor_lambda_snp)[2]='estimate'
# eff.FST.cor_lambda_snp$chrregion = rownames(eff.FST.cor_lambda_snp)
# eff.FST.cor_lambda_snp <- eff.FST.cor_lambda_snp %>%
#   mutate( confidence_interval = case_when(upper<0~"<0",
#                                           lower>0~">0",
#                                           upper>0 & lower<0 ~ "include 0")) %>%
#   mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0")))
# 
# 
# cor_lambda_fst_snp <- cor_lambda_fst_snp %>% inner_join(eff.FST.cor_lambda_snp %>% dplyr::select(chrregion, confidence_interval), by="chrregion")

# 
# 
# lalc <- lalc %>% 
#   filter(P1 != P2) %>%
#   mutate(chr=substr(chrregion,1,2)) %>%
#   mutate(region=substr(chrregion,3,5)) %>%
#   rename(correlation=estimateur) %>%
#   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
#   filter(region !="C") %>%
#   mutate(p=paste0(P1,"-",P2)) %>% 
#   mutate(p = factor(p, levels=c("WE-EE","EE-WA","WA-EA","WE-WA","EE-EA","WE-EA"))) %>%
#   droplevels() %>%
#   na.omit() %>%
#   dplyr::select(chrregion, chr, region, P1, P2, p, correlation, nb_intervalles, BIC) 
# 
# 
# dat1 <-  lalc %>%
#   inner_join(FSTb, by=c("chrregion","P1","P2")) %>%
#   na.omit() %>%
#   rename(correlation_common_SNP=correlation) %>%
#   inner_join(cor_lambda_fst_blocs, by=c("chrregion","P1","P2","FST","p","chr","region", "chrregion")) %>%
#   rename(correlation_specific_SNP=correlation) %>%
#   dplyr::select(P1,P2,FST,correlation_specific_SNP, correlation_common_SNP, chr, region, chrregion) 
# 
# 
# 
# 
# # SNP common versus specific
# mod = lmList( correlation_common_SNP ~ FST | chrregion, data=dat1)
# eff.FST.cor_lambda_fst_blocs_common_snp=data.frame(intervals(mod)[,,2])
# colnames(eff.FST.cor_lambda_fst_blocs_common_snp)[2]='estimate'
# eff.FST.cor_lambda_fst_blocs_common_snp$chrregion = eff.FST.cor_lambda_fst_blocs$chrregion
# eff.FST.cor_lambda_fst_blocs_common_snp$SNP = "common"
# 
# 
# eff.FST.cor_lambda_fst_blocs_snp_common <- lalc %>%
#   inner_join(FSTb, by=c("chrregion","P1","P2")) %>%
#   na.omit() %>%
#   inner_join(eff.FST.cor_lambda_fst_blocs_common_snp, by="chrregion") %>%
#   dplyr::select(chrregion, chr, region, p, correlation, FST,estimate, lower, upper)%>%
#   mutate( confidence_interval = case_when(upper<0~"<0",
#                                           lower>0~">0",
#                                           upper>0 & lower<0 ~ "include 0")) %>%
#   mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0")))
#   
# 
# length(which(eff.FST.cor_lambda_fst_blocs_common_snp$estimate <0))
# 
# eff.FST.cor_lambda_fst_blocs_common_snp %>%
#   mutate( confidence_interval = case_when(upper<0~"<0",
#                                           lower>0~">0",
#                                           upper>0 & lower<0 ~ "include 0")) %>%
#   mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0"))) %>%
#   dplyr::select(chrregion, confidence_interval) %>%
#   group_by(confidence_interval) %>%
#   summarise(n=n())
# 
# # check if there is a significant differences
# mod = lmList( correlation_common_SNP ~ FST | chrregion, data=dat1)
# common <- data.frame(summary(mod)$coefficients)[,c("Estimate.FST","Std..Error.FST")] %>%
#   rownames_to_column("chrregion") %>% 
#   rename(effect=Estimate.FST, sd=Std..Error.FST) %>% mutate(method="common") %>%
#   rowwise()
# 
# 
# 
# mod = lmList( correlation ~ FST | chrregion, data=cor_lambda_fst_blocs)
# specific <- data.frame(summary(mod)$coefficients)[,c("Estimate.FST","Std..Error.FST")] %>%
#   rownames_to_column("chrregion") %>%
#   rename(effect=Estimate.FST, sd=Std..Error.FST) %>% mutate(method="specific") %>%
#   rowwise()
# 
# tab <- rbind(common, specific) %>% mutate(chrregion=factor(chrregion), method=factor(method)) %>% mutate(c=1)
# 
# summary(lm(effect~method  ,weights=1/sd, data=tab))

# 
# hr <- fread(titre_hrec) %>% filter(kept==T)
# 
# 
# 
# 
# # Figure 4 : relationship correlation rho and CsRe at 4Mb in centromeres C
# 
# ratio=0.03
# # correlation rho-CsRe vs value of CsRe
# sf4ca <- carte_csre_landraces_4Mb %>% filter(region=="C") %>%
#   group_by(chr, region) %>%
#   summarise(correlation=cor(lambda_rho, cbay), mean_cbay=mean(cbay)) %>%
#   ggplot(aes(x=mean_cbay, y=correlation, col=region)) + geom_point() +
#   theme_light() +
#   geom_smooth(method = "lm", se=F) +
#   scale_colour_manual(values=C) +
#   ylab(expression(paste("correlation ", rho, " and CsRe")))  +
#   xlab("average CsRe Bayesian rec rates \n (cM/Mb)") +
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white")) +
#   coord_fixed(ratio=ratio)
# 
# 
# # correlation rho-CsRe vs value of rho
# ratio=0.0007
# sf4cb <- carte_csre_landraces_4Mb %>% filter(region=="C") %>%
#   group_by(chr, region) %>%
#   summarise(correlation=cor(lambda_rho, cbay), mean_lambda_rho=mean(lambda_rho)) %>%
#   ggplot(aes(x=mean_lambda_rho, y=correlation, col=region)) + geom_point() +
#   theme_light() +
#   geom_smooth(method = "lm", se=F) +
#   scale_colour_manual(values=C) +
#   ylab(expression(paste("correlation ", rho, " and CsRe")))  +
#   xlab("average LD-based rec rates \n (/kb)")  +
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white"))+
#   coord_fixed(ratio=ratio)
# 
# 
# 
# titre_sf4c <- paste0(dossier_graphes, "FigureS4C", extension)
# sf4c <- ggarrange(sf4ca,sf4cb, ncol=2, common.legend = T, legend="right", labels = c("A","B"),font.label = list(size = 18))
# sf4c
# 
# 
# tiff(titre_sf4c, units="in", width = 9, height = 9, res=200, compression = "lzw")
# sf4c
# dev.off()
# 
# 
# # Figure S4D : introgressions in CsRe
# # profiles in 7DR3 for CsRe
# sf4da <- carte_csre_landraces_4Mb %>% filter((chr=="7D" & region=="R3")) %>%
#   mutate(introgression=ifelse(posl >= 611e6, T, F)) %>%
#   ggplot(aes(x=posl, y=cbay, col=introgression, group=1)) + geom_line() +
#   ylab("CsRe Bayesian rec rates \n (cM/Mb)") +
#   xlab("physical position (Mb)") +
#   theme_light() +
#   ggtitle("7DR3")+
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white"))+
#   scale_x_continuous(labels=function(x)x/1e6)
# 
# # profiles in 7DR3 for WE
# sf4db <- carte_csre_landraces_4Mb %>% filter((chr=="7D" & region=="R3" & population=="WE")) %>%
#   mutate(introgression=ifelse(posl >= 611e6, T, F)) %>%
#   ggplot(aes(x=posl, y=lambda_rho, col=introgression , group=1)) + geom_line() +
#   ylab(expression(paste("WE ", rho, " (/kb)")))  +
#   xlab("physical position (Mb)") +
#   theme_light() +
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white"))+
#   scale_x_continuous(labels=function(x)x/1e6)
# 
# # scatterplot in 7DR3 CsRe vs WE rho
# sf4dc <- carte_csre_landraces_4Mb %>% filter((chr=="7D" & region=="R3" & population=="WE")) %>%
#   mutate(introgression=ifelse(posl >= 611e6, T, F)) %>%
#   ggplot(aes(x=cbay, y=lambda_rho, col=introgression , group=1)) + geom_point() +
#   ylab(expression(paste("WE ", rho, " (/kb)")))  +
#   xlab("CsRe Bayesian rec rates (cM/Mb)") +
#   theme_light() +
#   theme(legend.position = "none") +
#   geom_text(data=NULL, aes(x=Inf, y=0,  
#                            label=paste("cor = ", carte_csre_landraces_4Mb %>% filter((chr=="7D" & region=="R3" & population=="WE"))%>%
#                                          mutate(introgression=ifelse(posl >=611e6, T, F)) %>%
#                                          summarise(round(cor(lambda_rho, cbay), digits = 2)) %>%
#                                          unlist() %>% as.vector())),
#             size=5,
#             vjust=0,
#             hjust=+1.5,
#             col="black",
#             check_overlap = TRUE)+
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white"))
# 
# # profiles in 2AR1 for CsRe
# 
# sf4dd <- carte_csre_landraces_4Mb %>% filter((chr=="2A" & region=="R1")) %>%
#   mutate(introgression=ifelse(posl <=20e6, T, F)) %>%
#   ggplot(aes(x=posl, y=cbay, col=introgression, group=1)) + geom_line() +
#   ylab("CsRe Bayesian rec rates \n (cM/Mb)") +
#   xlab("physical position (Mb)") +
#   theme_light() +
#   ggtitle("2AR1") +
#   theme(plot.title = element_text(hjust = 0.5))+
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white"))+
#   scale_x_continuous(labels=function(x)x/1e6)
# 
# # profiles in 7DR3 for WE
# 
# sf4de <- carte_csre_landraces_4Mb %>% filter((chr=="2A" & region=="R1" & population=="WE")) %>%
#   mutate(introgression=ifelse(posl <=20e6, T, F)) %>%
#   ggplot(aes(x=posl, y=lambda_rho, col=introgression , group=1)) + geom_line() +
#   ylab(expression(paste("WE ", rho, " (/kb)")))  +
#   xlab("physical position (Mb)") +
#   theme_light() +
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white"))+
#   scale_x_continuous(labels=function(x)x/1e6)
# 
# 
# 
# # scatterplot in 2AR1 CsRe versus rho WE
# 
# sf4df <- carte_csre_landraces_4Mb %>% filter((chr=="2A" & region=="R1" & population=="WE")) %>%
#   mutate(introgression=ifelse(posl <=20e6, T, F)) %>%
#   ggplot(aes(x=cbay, y=lambda_rho, col=introgression , group=1)) + geom_point() +
#   ylab(expression(paste("WE ", rho, " (/kb)")))  +
#   xlab("CsRe Bayesian rec rates (cM/Mb)") +
#   theme_light() +
#   theme(legend.position = "none") +
#   geom_text(data=NULL, aes(x=Inf, y=0,  
#                            label=paste("cor = ", carte_csre_landraces_4Mb %>% filter((chr=="2A" & region=="R1" & population=="WE"))%>%
#                                          mutate(introgression=ifelse(posl <=20e6, T, F)) %>%
#                                          summarise(round(cor(lambda_rho, cbay), digits = 2)) %>%
#                                          unlist() %>% as.vector())),
#             size=5,
#             vjust=0,
#             hjust=+1.5,
#             col="black",
#             check_overlap = TRUE)+
#   theme(plot.title = element_text(hjust = 0.5, size=16),
#         strip.text.y = element_text(size = 14, color="black"),
#         axis.text.y = element_text(size=14),
#         axis.title.y = element_text(size=16),
#         axis.text.x = element_text(size=14),
#         axis.title.x = element_text(size=16),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=16),
#         strip.background = element_rect(color="black", fill="white"))
# 
# 
# sf4d <- ggarrange(sf4da, sf4dd, sf4db, sf4de, sf4dc, sf4df, ncol=2, nrow=3, common.legend = T, legend = "bottom", labels=c("A","D","B","E","C","F"),font.label = list(size = 18))
# sf4d
# 
# 
# 
# 
# tiff(titre_introgressions, units="in", width = 9, height = 12, res=200, compression = "lzw")
# sf4d
# dev.off()
# 
# 


# Figure 5 : boxplots of correlations for each pairs of pop

# significant differences in correlation







matrice_fst_blocs <- read.table(titre_matrice_fst_blocs) %>% as.matrix() %>%
  melt(value.name="FST") %>%
  rename(P1=Var1,P2=Var2) %>%
  mutate(P=paste0(P1,"-",P2)) %>%
  mutate(P = factor(P, levels=c("WE-EE","EE-WA","WA-EA","WE-WA","EE-EA","WE-EA"))) %>%
  na.omit() %>%
  mutate(FST=100*round(FST,3)) %>%
  arrange(FST) %>%
  dplyr::select(P, FST)

# if (repertoire !="common_SNP"){

fread(titre_correlations_mixed_models_published) %>% 
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
  arrange(signif)%>%
  dplyr::select(P, signif)



#}



f4 <- fread(titre_correlations_mixed_models_published) %>% 
  mutate(P=paste0(P1,"-",P2)) %>%
  mutate(P = factor(P, levels=c("WE-EE","WE-WA", "EE-WA","WA-EA","WE-EA", "EE-EA"))) %>%
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
tiff(titre_f4, units="in", width = 12, height = 9, res=200, compression = "lzw")
f4f
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))+
  theme(legend.position="bottom") +
  ylab(expression(paste("correlation of log10(", lambda,")"))) +
  xlab("Fst (%, from haplotypic alleles)") +
  guides(colour=guide_legend("95% confidence interval of slope estimate : "))+
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+
  scale_x_continuous( labels=function(x)x*100,breaks=seq(0,1,0.05))+
  theme(strip.text.y = element_text(size = 14, color="black"),
        strip.text.x = element_text(size = 14, color="black"),
        strip.background = element_rect(color="black", fill="white"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
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
  theme( axis.text.y = element_text(size=6)) +
  theme(legend.position="none") +
  guides(colour=guide_legend("95% confidence interval of slope estimate: "))+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        panel.grid.minor = element_blank())




f6 <- ggdraw() +
  draw_plot(f6a, x=-0.33, y=0,1, 0.5, width =1.35) +
  draw_plot(f6b,x=0.7,y=0,1,width=0.3)+
  draw_plot_label(c("A", "B"), c(0,0.68), c(0.8,1), size = 20)


f6


tiff(titre_f6, units="in", width = 18, height = 9, res=200, compression = "lzw")
f6
dev.off()





# 
# 
# 
# 
# # relationship correlation lambda fst in 3BR1
# example_correlation_lambda_fst_3BR1 <- cor_lambda_fst_blocs %>%
#   filter(chr=="3B" & region=="R1") %>%
#   ggplot(aes(x=FST, y=correlation, col=p)) + 
#   geom_point(size=5) +
#   geom_smooth(method = "lm",se=F, col="#F8766D", size=2) +
#   theme_light()+
#   scale_colour_manual(values=c(viridis(t)[n+6],viridis(t)[n+5],viridis(t)[n+4],viridis(t)[n+3],viridis(t)[n+2],viridis(t)[n+1]))+
#   scale_x_continuous( labels=function(x)x*100)+
#   theme(strip.text.y = element_text(size = 18, color="black"),
#         strip.text.x = element_text(size = 18, color="black"),
#         strip.background = element_rect(color="black", fill="white"),
#         axis.title.y = element_text(size=18),
#         axis.title.x = element_text(size=18),
#         axis.text.y = element_text(size=16),
#         axis.text.x = element_text(size=16),
#         legend.text=element_text(size=18),
#         legend.title=element_text(size=0),
#         legend.key.size = unit(1, "cm"),
#         plot.title = element_text(hjust = 0.5, size=20),
#         panel.grid.minor = element_blank())+
#   ylab(expression(paste("correlation of log10(", lambda,")"))) +
#   xlab("Fst (%, from haplotypic alleles)") +
#   ggtitle("3BR1")
# 
# tiff(titre_example_correlation_lambda_fst_3BR1, units="in", width = 12, height = 9, res=200, compression = "lzw")
# example_correlation_lambda_fst_3BR1
# dev.off()
# 
# 
# f6ea <- cor_lambda_fst_snp %>%
#   ggplot(aes(x=FST, y=correlation)) + 
#   geom_point() +
#   geom_smooth(method="lm", se=F, aes(col=confidence_interval)) +
#   scale_color_manual(values=c("#F8766D","#619CFF","#00BA38")) +
#   
#   facet_grid(region~chr, scales="free") +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))+
#   theme(legend.position="bottom") +
#   ylab(expression(paste("correlation of log10(", lambda,")"))) +
#   xlab("Fst (%, from SNP alleles)") +
#   guides(colour=guide_legend("95% confidence interval of slope estimate : "))+
#   scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+
#   scale_x_continuous( labels=function(x)x*100,breaks=seq(0,1,0.1))+
#   theme(strip.text.y = element_text(size = 14, color="black"),
#         strip.text.x = element_text(size = 14, color="black"),
#         strip.background = element_rect(color="black", fill="white"),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=14),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=14),
#         panel.grid.minor = element_blank())
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# f6eb <- eff.FST.cor_lambda_snp %>%
#   mutate( cr = fct_reorder(chrregion, estimate)) %>%
#   ggplot(aes(x=estimate,y=cr,colour=confidence_interval)) +
#   scale_color_manual(values=c("#F8766D","#619CFF","#00BA38")) +
#   
#   geom_segment(aes(x=lower,xend=upper,yend=chrregion),colour='gray', size=1) +
#   geom_point() +
#   geom_vline(xintercept=0,col='#619CFF', size=1) +
#   xlab(label='Fst effect') + ylab( label = "genomic region") +
#   theme_light()+
#   theme( axis.text.y = element_text(size=6)) +
#   theme(legend.position="none") +
#   guides(colour=guide_legend("95% confidence interval of slope estimate: "))+
#   theme(axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=14),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=14),
#         panel.grid.minor = element_blank())
# 
# 
# 
# f6e <- ggdraw() +
#   draw_plot(f6ea, x=-0.33, y=0,1, 0.5, width =1.35) +
#   draw_plot(f6eb,x=0.7,y=0,1,width=0.3)+
#   draw_plot_label(c("A", "B"), c(0,0.68), c(0.8,1), size = 20)
# 
# 
# f6e
# 
# 
# tiff(titre_f6e, units="in", width = 18, height = 9, res=200, compression = "lzw")
# f6e
# dev.off()
# 
# 
# 
# 
# 
# 
# # # relation lambda fst with common SNP
# # 
# # f6acommon <- eff.FST.cor_lambda_fst_blocs_snp_common %>% 
# #   ggplot(aes(x=FST, y=correlation)) + 
# #   geom_point() +
# #   geom_smooth(method="lm", se=F, aes(col=confidence_interval)) +
# #   scale_color_manual(values=c("#F8766D","#619CFF","#00BA38")) +
# #   facet_grid(region~chr, scales="free") +
# #   theme_light() +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))+
# #   theme(legend.position="bottom") +
# #   ylab(expression(paste("correlation of log10(", lambda,")"))) +
# #   xlab("Fst (%, from SNP alleles)") +
# #   guides(colour=guide_legend("95% confidence interval of slope estimate : "))+
# #   scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+
# #   scale_x_continuous( labels=function(x)x*100,breaks=seq(0,1,0.1))+
# #   theme(strip.text.y = element_text(size = 14, color="black"),
# #         strip.text.x = element_text(size = 14, color="black"),
# #         strip.background = element_rect(color="black", fill="white"),
# #         axis.title.y = element_text(size=14),
# #         axis.title.x = element_text(size=14),
# #         legend.text=element_text(size=14),
# #         legend.title=element_text(size=14),
# #         panel.grid.minor = element_blank())
# #   theme(plot.title = element_text(hjust = 0.5)) 
# #   
# # f6bcommon <- eff.FST.cor_lambda_fst_blocs_snp_common %>%
# #     mutate( cr = fct_reorder(chrregion, estimate)) %>%
# #     ggplot(aes(x=estimate,y=cr,colour=confidence_interval)) +
# #     scale_color_manual(values=c("#F8766D","#619CFF","#00BA38")) +
# #     geom_segment(aes(x=lower,xend=upper,yend=chrregion),colour='gray', size=1) +
# #     geom_point() +
# #     geom_vline(xintercept=0,col='#619CFF', size=1) +
# #     xlab(label='Fst effect') + ylab( label = "genomic region") +
# #     theme_light()+
# #     theme( axis.text.y = element_text(size=6)) +
# #     theme(legend.position="none") +
# #     guides(colour=guide_legend("95% confidence interval of slope estimate: "))+
# #     theme(axis.title.y = element_text(size=14),
# #           axis.title.x = element_text(size=14),
# #           legend.text=element_text(size=14),
# #           legend.title=element_text(size=14),
# #           panel.grid.minor = element_blank())
# 
# # 
# # 
# # sf6da <- dat1 %>% 
# #   ggplot(aes(x=FST, y=correlation_common_SNP)) + 
# #   geom_point() +
# #   geom_smooth(method="lm", se=F, col="#F1A340") +
# #   facet_grid(region~chr, scales = "free") +
# #   labs(colour = NULL) +
# #   theme_light()  +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))+
# #   theme(legend.position="none") +
# #   ylab(expression(paste("correlation of log10(", lambda,")"))) +
# #   xlab("Fst (%, from haplotypic alleles)") +
# #   scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+
# #   scale_x_continuous( labels=function(x)x*100,breaks=c(0,0.05,0.1,0.15, 0.2))+
# #   theme(strip.text.y = element_text(size = 14, color="black"),
# #         strip.text.x = element_text(size = 14, color="black"),
# #         strip.background = element_rect(color="black", fill="white"),
# #         axis.title.y = element_text(size=14),
# #         axis.title.x = element_text(size=14),
# #         legend.text=element_text(size=14),
# #         legend.title=element_text(size=14),
# #         panel.grid.minor = element_blank())+
# #   ggtitle("Common SNP dataset") + 
# #   theme(plot.title = element_text(hjust = 0.5))
# # 
# # 
# # 
# # 
# # sf6db <-  dat1 %>% 
# #   ggplot(aes(x=FST, y=correlation_specific_SNP)) + 
# #   geom_point() +
# #   geom_smooth(method="lm", se=F, col="#C77CFF") +
# #   facet_grid(region~chr, scales = "free") +
# #   labs(colour = NULL) +
# #   theme_light() +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))+
# #   theme(legend.position="none") +
# #   ylab(expression(paste("correlation of log10(", lambda,")"))) +
# #   xlab("Fst (%, from haplotypic alleles)") +
# #   scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+
# #   scale_x_continuous( labels=function(x)x*100,breaks=c(0,0.05,0.1,0.15, 0.2))+
# #   theme(strip.text.y = element_text(size = 14, color="black"),
# #         strip.text.x = element_text(size = 14, color="black"),
# #         strip.background = element_rect(color="black", fill="white"),
# #         axis.title.y = element_text(size=14),
# #         axis.title.x = element_text(size=14),
# #         legend.text=element_text(size=14),
# #         legend.title=element_text(size=14),
# #         panel.grid.minor = element_blank())+
# #   ggtitle("Population specific SNP dataset") + 
# #   theme(plot.title = element_text(hjust = 0.5))
# # 
# # 
# # 
# # sf6dc <- eff.FST.cor_lambda_fst_blocs  %>%
# #   mutate(SNP="specific") %>%
# #   dplyr::select(-confidence_interval) %>%
# #   rbind(., eff.FST.cor_lambda_fst_blocs_common_snp) %>%
# #   mutate(SNP=factor(SNP, levels=c("specific","common"))) %>%
# #   group_by(chrregion) %>% 
# #   mutate(m=mean(estimate)) %>%
# #   ungroup() %>%
# #   mutate( cr = fct_reorder(chrregion, m)) %>%
# #   ggplot(aes(y = estimate, ymin = lower,
# #              x = cr, ymax = upper, group = row.names(.), col=SNP)) +
# #   geom_point(position = position_dodge(.5), size=4) +
# #   geom_linerange(position = position_dodge(.5)) +
# #   coord_flip()   +
# #   geom_hline(yintercept=0,col='#00BFC4', size=1) +
# #   ylab(label='FST effect') + xlab( label = "genomic region") +
# #   theme_minimal()+
# #   scale_colour_manual(values=c("#C77CFF","#F1A340"))+
# #   theme(plot.margin = unit(c(1,0,0,0), "cm")) +
# #   theme_light()+
# #   theme( axis.text.y = element_text(size=6)) +
# #   guides(colour=guide_legend("SNP dataset"))+
# #   theme(axis.title.y = element_text(size=18),
# #         axis.text.x = element_text(size=18),
# #         axis.title.x = element_text(size=18),
# #         legend.text=element_text(size=18),
# #         legend.title=element_text(size=18),
# #         panel.grid.minor = element_blank(),
# #         legend.position = "bottom")
# # 
# # sf6dd <- eff.FST.cor_lambda_fst_blocs  %>%
# #   mutate(SNP="specific") %>%
# #   dplyr::select(-confidence_interval) %>%
# #   rbind(., eff.FST.cor_lambda_fst_blocs_common_snp) %>%
# #   pivot_wider(id_cols = "chrregion", "SNP", values_from = c("lower","estimate","upper")) %>%
# #   mutate(signif=case_when(upper_specific <0 & upper_common < 0 ~ "< 0 in both analysis",
# #                           (upper_specific >0 & upper_common < 0) |  (upper_specific <0 & upper_common > 0)~ "include 0 in at least one analysis",
# #                           (upper_specific >0 & upper_common < 0) |  (upper_specific <0 & upper_common > 0)~ "include 0 in at least one analysis",
# #                           lower_specific < 0 & lower_common < 0 & (upper_specific >0 | upper_common > 0)~ "include 0 in at least one analysis")) %>%
# #     ggplot(aes(x=estimate_specific,y=estimate_common, col=signif)) + geom_point(size=4) +
# #   guides(colour=guide_legend("95% confidence interval of slopes estimates: "))+
# #   theme_light()+
# #   
# #   theme(axis.title.y = element_text(size=18),
# #         axis.title.x = element_text(size=18),
# #         axis.text.x = element_text(size=18),
# #         axis.text.y = element_text(size=18),
# #         
# #         legend.text=element_text(size=18),
# #         legend.title=element_text(size=18),
# #         panel.grid.minor = element_blank(),
# #         legend.position = "bottom") +
# #   xlab("Slope estimate with specific SNP dataset") +
# #   ylab("Slope estimate with common SNP dataset") +
# #   geom_abline(slope=1, intercept = 0, col="black") 
# #   
# # # 
# # # g <- ggdraw() +
# # #   draw_plot(g2, x=0, y=0.5, width = 0.5, height = 0.5)+
# # #   draw_plot(g1, x=0.5, y=0.5, width = 0.5, height = 0.5)+
# # #   draw_plot(g3, x=0, y=0, width = 1, height = 0.5) +
# # #   draw_plot_label(c("A","B","C"), c(0,0.5,0), c(1, 1, 0.5), size = 15)
# # 
# # 
# # 
# # sf6d <- ggdraw() +
# #   draw_plot(sf6dc, x=0, y=0.5, width = 1, height = 0.5)+
# #   draw_plot(sf6dd, x=0, 0, width =1, height = 0.5)+
# #   draw_plot_label(c("A","B"), c(0,0), c(1, 0.5), size = 15)
# # 
# # sf6d
# # 
# # 
# # tiff(titre_sf6d, units="in", width = 14, height = 16, res=200, compression = "lzw")
# # sf6d
# # dev.off()
# 
# 
# 
# 
# # Relationship correlation of rho and fst
# rho <- fread(titre_rho_aleatoires) %>% 
#   filter(P1 != P2) %>%
#   mutate(chr=substr(chrregion,1,2)) %>%
#   mutate(region=substr(chrregion,3,5)) %>%
#   rename(correlation=estimateur) %>%
#   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
#   filter(region !="C") %>%
#   mutate(p=paste0(P1,"-",P2)) %>% 
#   mutate(p = factor(p, levels=c("WE-EE","EE-WA","WA-EA","WE-WA","EE-EA","WE-EA"))) %>%
#   droplevels() %>%
#   na.omit() %>%
#   dplyr::select(chrregion, chr, region, P1, P2, p, correlation, nb_intervalles, BIC) 
# 
# 
# 
# dat <-  rho %>%
#   inner_join(FSTb, by=c("chrregion","P1","P2")) %>%
#   na.omit()
# 
# 
# 
# mod = lmList( correlation ~ FST | chrregion, data=dat)
# eff.FST.cor_rho_fst_haplotypic=data.frame(intervals(mod)[,,2])
# colnames(eff.FST.cor_rho_fst_haplotypic)[2]='estimate'
# eff.FST.cor_rho_fst_haplotypic$chrregion = rownames(eff.FST.cor_rho_fst_haplotypic)
# eff.FST.cor_rho_fst_haplotypic <- eff.FST.cor_rho_fst_haplotypic %>%
#   mutate( confidence_interval = case_when(upper<0~"<0",
#                                           lower>0~">0",
#                                           upper>0 & lower<0 ~ "include 0")) %>%
#   mutate(confidence_interval=factor(confidence_interval, levels=c("<0","include 0",">0")))
# 
# 
# dat <- dat %>% inner_join(eff.FST.cor_rho_fst_haplotypic %>% dplyr::select(chrregion, confidence_interval), by="chrregion")
# 
# g1 <-dat %>%
#   ggplot(aes(x=FST, y=correlation)) + 
#   geom_point() +
#   geom_smooth(method="lm", se=F, aes(col=confidence_interval)) +
#   facet_grid(region~chr, scales = "free") +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))+
#   theme(legend.position="bottom") +
#   ylab(expression(paste("correlation of log10(", rho,")"))) +
#   xlab("Fst (%, from haplotypic alleles)") +
#   guides(colour=guide_legend("95% confidence interval of slope estimate : "))+
#   scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1))+
#   scale_x_continuous( labels=function(x)x*100,breaks=c(0,0.05,0.1,0.15, 0.2))+
#   theme(strip.text.y = element_text(size = 14, color="black"),
#         strip.text.x = element_text(size = 14, color="black"),
#         strip.background = element_rect(color="black", fill="white"),
#         axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=14),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=14),
#         panel.grid.minor = element_blank())+
#   scale_colour_manual(values=c("#F8766D","#00BFC4", "#00BA38"))
# 
# 
# 
# 
# g2 <- eff.FST.cor_rho_fst_haplotypic %>%
#   mutate( cr = fct_reorder(chrregion, estimate)) %>%
#   ggplot(aes(x=estimate,y=cr,colour=confidence_interval)) +
#   geom_segment(aes(x=lower,xend=upper,yend=chrregion),colour='gray') +
#   geom_point() +
#   geom_vline(xintercept=0,col='#00BFC4', size=1) +
#   xlab(label='Fst effect') + ylab( label = "genomic region") +
#   theme_light()+
#   theme( axis.text.y = element_text(size=6)) +
#   theme(legend.position="none") +
#   guides(colour=guide_legend("95% confidence interval of slope estimate: "))+
#   theme(axis.title.y = element_text(size=14),
#         axis.title.x = element_text(size=14),
#         legend.text=element_text(size=14),
#         legend.title=element_text(size=14),
#         panel.grid.minor = element_blank())+
#   scale_colour_manual(values=c("#F8766D","#00BFC4", "#00BA38"))
# 
# 
# g <- ggarrange(g1,g2,ncol=2, common.legend = T, labels=c("A","B"), legend = "bottom", widths = c(2,1))
# 
# table(eff.FST.cor$confidence_interval)
# eff.FST.cor[which(eff.FST.cor$estimate>0),]
# eff.FST.cor[which(eff.FST.cor$estimate<0),] %>%
#   summarise(min=min(estimate), max=max(estimate))
# 
# g
# ggsave(titre_correlation_rho_FST_haplotypic_blocks,g, width = 40, height = 20, units="cm")
# 
# 
# 
# 
# 
# ### 3) Pairwise correlations of lambda
# 
# # table(lal$BIC)/6
# 
# lal %>% summarise(nb_intervalles_min=min(nb_intervalles), nb_intervalles_max=max(nb_intervalles))
# lal %>% dplyr::select(chrregion, BIC) %>%
#   unique() %>%
#   group_by(BIC) %>% 
#   summarise(n=n())
# 
# pairwise.t.test(lal$correlation, lal$p, p.adjust.method = "bonf")
# 
# 
# 
# 
# lal %>% summarise(correlation_lambda_moyen=round(mean(correlation), digits=2), sd=round(sd(correlation),2 )) 
# 
# 
# 
# # g <- ggplot(lal, aes(x=p, y=correlation, col=p)) + geom_boxplot()  +
# #   theme_light() +
# #   theme(legend.position = "none")+
# #   annotate("text", x = "WE-EE", y = 1.1, label = "a") +
# #   annotate("text", x = "WE-WA", y = 1.1, label = "b") +
# #   annotate("text", x = "EE-WA", y = 1.1, label = "b") +
# #   annotate("text", x = "WA-EA", y = 1.1, label = "c") +
# #   annotate("text", x = "WE-EA", y = 1.1, label = "cd") +
# #   annotate("text", x = "EE-EA", y = 1.1, label = "d") +
# #   xlab("") +
# #   ylab(expression(paste("Correlation log10(", lambda, ")")))
# 
# 
# # texte <- cbind(data.frame(P1=combn(c("WE","EE","WA","EA"), m=2)),data.frame(P1=combn(rev(c("WE","EE","WA","EA")), m=2))) %>%
# #   t() %>%
# #   as.data.frame() %>%
# #   rename(P1=V1,P2=V2) %>%
# #   filter(P1 !=P2) %>%
# #   mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
# #   mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
# #   arrange(P1,P2) %>%
# #   mutate(signif=c("a","b","cd","a","b","d","b","b","c","cd","d","c"))
# 
# 
# 
# 
# #WE-EE : a
# #WE-WA : ab
# #EE-WA : ab
# #WA-EA : bc
# #WE-EA : c
# #EE-EA : c
# 
# # texte <- cbind(data.frame(P1=combn(c("WE","EE","WA","EA"), m=2)),data.frame(P1=combn(rev(c("WE","EE","WA","EA")), m=2))) %>%
# #   t() %>%
# #   as.data.frame() %>%
# #   rename(P1=V1,P2=V2) %>%
# #   filter(P1 !=P2) %>%
# #   mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
# #   mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
# #   arrange(P1,P2) %>%
# #   mutate(signif=c("a","b","cd","a","b","d","b","b","c","cd","d","c"))
# 
# 
# 
# 
# # 
# # n=4
# # t <- 10
# # g <- ggplot(lal2, aes(x=factor(1), y=correlation)) + 
# #   geom_boxplot() +
# #   facet_grid(P1~P2) +
# #   theme_light() +
# #   theme(axis.title.x=element_blank(),
# #         axis.text.x=element_blank(),
# #         axis.ticks.x=element_blank(),
# #         legend.position = "none") +
# #   geom_text(
# #     data    = texte,
# #     mapping = aes(x = -Inf, y = Inf, label = signif),
# #     hjust   = -1,
# #     vjust   = +1.5
# #   ) +
# #   ylab(expression(paste("Correlation log10(", lambda, ")")))+
# #   geom_boxplot(data=lal2, mapping=aes(x=factor(1), y=correlation, fill=color, middle = mean(correlation))) +
# #   scale_fill_manual(values=c(viridis(t)[n+3],viridis(t)[n+2],viridis(t)[n+4],viridis(t)[n+1],viridis(t)[n+5],viridis(t)[n+6])) 
# # 
# # titre_g <- paste0(r_graphe,"lambda_pairwise", extension)
# # 
# # g
# # ggsave(titre_g, g)
# 
# 
# 
# 
# 
# 
# #### fst compute on SNP
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Correlation of 4 Mb windows genomewide per population
# carte %>% group_by(population) %>%
#   summarise(correlation=round(cor(cbay, lambda_rho), digits=2))
# 
# # Significant differences in correlation coefficients between populations
# data <- carte %>% group_by(population) %>%
#   pivot_wider(id_cols = c("chr","posl","posr","cbay"), names_from = "population", values_from = "lambda_rho") %>%
#   na.omit() %>%
#   mutate(cbay=as.numeric(cbay)) %>%
#   as.data.frame()
# cocor(~ WE +cbay | EE + cbay, data=data, return.htest = T, test="zou2007")
# cocor(~ WE +cbay | WA + cbay, data=data, return.htest = T, test="zou2007")
# cocor(~ WE +cbay | EA + cbay, data=data, return.htest = T, test="zou2007")
# cocor(~ EE +cbay | WA + cbay, data=data, return.htest = T, test="zou2007")
# cocor(~ EE +cbay | EA + cbay, data=data, return.htest = T, test="zou2007")
# cocor(~ WA +cbay | EA + cbay, data=data, return.htest = T, test="zou2007")
# 
# 
# 
# # Correlation of 4 Mb windows averaged over all populations and regions
# carte %>% group_by(population, chr, region) %>%
#   summarise(correlation=cor(lambda_rho, cbay), nb_intervalles=n()) %>%
#   ungroup() %>%
#   filter(nb_intervalles >= 5) %>%
#   summarise(mean_correlation_genomic_region=round(mean(correlation, na.rm=T), digits=2), 
#             sd_correlation_genomic_region=round(sd(correlation, na.rm=T), digits=2))
# 
# # Correlation of 4 Mb windows averaged over all populations and for C region only
# carte %>% group_by(population, chr, region) %>%
#   filter(region=="C") %>%
#   summarise(correlation=cor(lambda_rho, cbay), nb_intervalles=n()) %>%
#   ungroup() %>%
#   filter(nb_intervalles >= 5) %>%
#   summarise(mean_correlation_genomic_region=round(mean(correlation, na.rm=T), digits=2), 
#             sd_correlation_genomic_region=round(sd(correlation, na.rm=T), digits=2))
# carte %>% group_by(population, chr, region) %>%
#   filter(region!="C") %>%
#   summarise(correlation=cor(lambda_rho, cbay), nb_intervalles=n()) %>%
#   ungroup() %>%
#   filter(nb_intervalles >= 5) %>%
#   summarise(mean_correlation_genomic_region=round(mean(correlation, na.rm=T), digits=2), 
#             sd_correlation_genomic_region=round(sd(correlation, na.rm=T), digits=2))
# 
# 
# 
# 
# test <- carte %>% group_by(population, chr, region ) %>%
#   summarise(correlation=cor(lambda_rho, cbay), nb_intervalles=n()) %>%
#   filter(nb_intervalles>=5) 
# 
# pairwise.wilcox.test(test$correlation, test$population, p.adjust.method = "bonf")
# pairwise.t.test(test$correlation, test$population, p.adjust.method = "bonf")
# 
# test %>% group_by(population) %>%
#   summarise(m=mean(correlation),s=sd(correlation))
# 
# 
# # centromeres
# 
# data <- carte %>% filter(region=="C") %>%
#   group_by(chr, region) %>%
#   summarise(correlation=cor(lambda_rho, cbay), mean_cbay=mean(cbay), mean_lambda_rho=mean(lambda_rho)) 
# summary(lm(correlation~mean_cbay, data=data))
# summary(lm(correlation~mean_lambda_rho, data=data))
# 
# 
# 
# ### Comparison of recombination rates between the most diverging populations at 4 Mb
# carte %>% pivot_wider(id_cols = c("chr","posl","posr", "region"), values_from = "lambda_rho", names_from = "population") %>%
#   na.omit() %>%
#   dplyr::select(chr, region, WE,EE,WA,EA) %>%
#   group_by(chr, region) %>%
#   summarise(WEEE=cor(WE,EE),
#             WEWA=cor(WE,WA),
#             WEEA=cor(WE,EA),
#             EEWA=cor(EE,WA),
#             EEEA=cor(EE,EA),
#             WAEA=cor(WA,EA)) %>%
#   pivot_longer(cols=c("WEEE","WEWA","WEEA","EEWA","EEEA","WAEA")) %>%
#   group_by(name) %>%
#   summarise(m=mean(value),s=sd(value)) %>%
#   arrange(desc(m))
# 
# #### 2) Introgressions
# 
# # Correlation of rho and CsRe per 4 Mb windows averaged over all population for regions of introgressions only
# carte %>% filter((region=="R3" & chr=="7D")|(region=="R1" & chr=="2A")) %>%
#   group_by(population, chr, region) %>%
#   summarise(correlation=cor(lambda_rho, cbay)) %>%
#   group_by(chr, region) %>%
#   summarise(mean_correlation=round(mean(correlation), digits=2),
#             sd_correlation=round(sd(correlation), digits=2))
# 
# 
# 
# 
# # variance of estimator
# 
# # 
# # carte <- csre_4Mb %>% 
# #   dplyr::select(chr, posl,posr,cbay, sdCbay) %>% 
# #   inner_join(landraces_4Mb, by=c("chr","posl","posr")) %>%
# #   dplyr::select(population, chr, region, posl, posr, cbay, lambda_rho, sdCbay, sd_lambda_rho) %>%
# #   arrange(population, chr, posl) %>%
# #   mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
# #   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
# #   mutate(lambda_rho=lambda_rho*1e3) %>%
# #   na.omit() %>%
# #   group_by(population, chr, region) %>%
# #   summarise(correlation=cor(cbay, lambda_rho), mean_sdCbay=mean(sdCbay), mean_sd_lambda_rho=mean(sd_lambda_rho)) %>%
# #   na.omit()
# # 
# # 
# # carte %>% filter(population=="WE" & region=="C") %>%ggplot(aes(x=mean_sdCbay, y=correlation, col=region)) +
# #   geom_point() + geom_smooth(method = "lm") +
# #   theme_light() +
# #   xlab("average standard error of CsRe rec rate on 4 Mb windows") +
# #   ylab("correlation")
# # 
# # 
# # carte %>% filter(population=="WE" & region=="C") %>%
# #   ggplot(aes(x=mean_sd_lambda_rho, y=correlation, col=region)) +
# #   geom_point() + geom_smooth(method = "lm") +
# #   theme_light() +
# #   xlab("average standard deviation of WE rec rate on 4 Mb windows") +
# #   ylab("correlation")
# # 
# # carte %>% filter(population=="WE") %>%
# #   ggplot(aes(x=mean_sd_lambda_rho, y=mean_sdCbay, col=region)) + geom_point() + geom_smooth(method = "lm") +
# #   theme_light() +
# #   xlab("average standard deviation of lambda_rho") +
# #   ylab("average standard error of Cbay") +
# #   scale_x_continuous(trans="log10") +
# #   scale_y_continuous(trans = "log10")
# 
# # correlation at 4Mb scale genomewide
# carte %>% pivot_wider(names_from = "population", values_from = "lambda_rho") %>%
#   na.omit() %>%
#   summarise(WEEA=round(cor(WE,EA, method = "spearman"), digits=2),
#             EEEA=round(cor(EE,EA, method = "spearman"), digits=2),
#             WEEE=round(cor(WE,EE, method = "spearman"), digits=2))
# 
# 
# 
# 
# 
# ### 3 ) Relationship lambda and Fst
# 
# 
# cat("\n\n Estimation of correlation of historical profiles within genomic regions of 3B chromosome computed on common SNP dataset \n")
# lalc <- fread(titre_cor_lambda_SNP_communs)
# 
# cat("\n\n Fst within genomic regions computed with haplotypic alleles \n")
# FSTb <- fread(titre_FST_blocs) 
# head(FSTb)
# 
# cat("\n\n Fst within genomic regions computed with SNP alleles \n")
# FSTs <- fread(titre_FST_SNP) 
# head(FSTs)
# 
# 
# 
# cor_lambda_fst_blocs <-  lal %>%
#   inner_join(FSTb, by=c("chrregion","P1","P2")) %>%
#   na.omit()
# # 
# # g <-cor_lambda_fst_blocs %>%
# #   ggplot(aes(x=FST, y=correlation)) + 
# #   geom_point() +
# #   geom_smooth(method="lm", se=F, col="black") +
# #   facet_grid(region~chr, scales = "free") +
# #   geom_point(data=cor_lambda_fst_blocs, aes(x=FST, y=correlation, col=p)) +
# #   labs(colour = NULL) +
# #   theme_light() +
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))+
# #   theme(legend.position="bottom") +
# #   guides(colour = guide_legend(nrow = 1))  +
# #   ylab(expression(paste("Correlation of log10(", lambda,")"))) +
# #   xlab("Fst (haplotypic alleles)") 
# # 
# # titre_graphe <- paste0(r_graphe, "lambda_Fst_haplotypes_populations", extension)
# # g
# # ggsave(titre_graphe, g)
# 
# 
# 
# 
# 
# 
# 
# 
# # FST compute with SNP
# 
# # Common SNP dataset
# 
# 
# # 
# # 
# # titre_g <- paste0(r_graphe,"common_or_specific_SNP_dataset", extension)
# # ggsave(titre_g,g, width = 40, height = 20, units="cm")
# 
# temp <- eff.FST.cor2  %>%
#   rbind(., eff.FST.cor) %>%
#   mutate(SNP=factor(SNP, levels=c("specific","common"))) %>%
#   pivot_wider(id_cols = chrregion, names_from = c("SNP"), values_from = c("lower","upper")) %>%
#   mutate(significant_differences=ifelse(!(upper_specific <= lower_common | upper_common <= lower_specific) ==T, FALSE,TRUE)) %>%
#   dplyr::select(chrregion, significant_differences)
# table(temp$significant_differences)
# temp[which(temp$significant_differences==T),]
# 
# 
# ### 4) SNP density
# 
# cat("\n\n summary outputs from PHASE \n")
# res <- fread(titre_resume) %>% filter(w_center==T)
# 
# 
# res %>% dplyr::select(population, SNPlpop, SNPrpop, lpop) %>%
#   unique() %>%
#   group_by(population) %>%
#   summarise(nb_SNP=n(), d_phy=sum(lpop), d_SNP=nb_SNP/d_phy) %>%
#   ungroup() %>%
#   summarise(d_SNP=mean(d_SNP))
# 
# 
# res %>% dplyr::select(population, region, SNPlpop, SNPrpop, lpop) %>%
#   unique() %>%
#   group_by(population, region) %>%
#   summarise(nb_SNP=n(), d_phy=sum(lpop), d_SNP=nb_SNP/d_phy) %>%
#   group_by(region) %>%
#   summarise(d_SNP=mean(d_SNP))
# 
# res %>% summarise(l=mean(lpop))
# 
# 
# res %>%
#   mutate(groupe=substr(chr,2,2)) %>%
#   dplyr::select(population, groupe, SNPlpop, SNPrpop, lpop) %>%
#   unique() %>%
#   group_by(population, groupe) %>%
#   summarise(nb_SNP=n(), d_phy=sum(lpop), d_SNP=nb_SNP/d_phy) %>%
#   group_by(groupe) %>%
#   summarise(d_SNP=mean(d_SNP))
# 
# 
# 
# 
# 
# ### 7) hrec
# 
# cat("\n\n HR intervals \n")
# hrec <- fread(titre_hrec) 
# head(hrec)
# 
# hrec2 <- hrec %>% filter(kept==T) %>%
#   dplyr::select(population, chr, region, IDintHR) %>% 
#   unique()
# 
# round(table(hrec2$population)/nrow(hrec2), digits=2)
# round(table(hrec2$region)/nrow(hrec2), digits=3)
# 
# table(hrec2$population, hrec2$chr, hrec2$region)
# 
# 
# fread(titre_resume) %>% filter() %>%
#   filter(w_center==T) %>%
#   filter(region !="C") %>%
#   na.omit() %>%
#   group_by(population) %>%
#   summarise(nbint=n()) %>%
#   cbind(as.data.frame(table(hrec2$population))) %>%
#   mutate(proportion_intervals_per_pop = round(Freq/nbint, digits=3))
# 
# 
# pop <- fread(titre_resume) %>% filter() %>%
#   filter(w_center==T) %>%
#   filter(region !="C") %>%
#   na.omit() %>%
#   group_by(population) %>%
#   summarise(nbint_pop=n()) %>%
#   cbind(as.data.frame(table(hrec2$region))) %>%
#   dplyr::select(nbint_pop, Freq) %>%
#   as.matrix()
# 
# pop
# test <- chisq.test(pop)
# test
# test$expected
# 
# 
# 
# fread(titre_resume) %>% filter() %>%
#   filter(w_center==T) %>%
#   filter(region !="C") %>%
#   na.omit() %>%
#   group_by(region) %>%
#   summarise(nbint=n()) %>%
#   cbind(as.data.frame(table(hrec2$region))) %>%
#   mutate(proportion_intervals_per_region = round(Freq/nbint, digits=3))
# 
# 
# 
# regions <- fread(titre_resume) %>% filter() %>%
#   filter(w_center==T) %>%
#   filter(region !="C") %>%
#   na.omit() %>%
#   group_by(region) %>%
#   summarise(nbint_region=n()) %>%
#   cbind(as.data.frame(table(hrec2$region))) %>%
#   dplyr::select(nbint_region, Freq) %>%
#   as.matrix()
# 
# regions
# test <- chisq.test(regions)
# test
# test$expected
# 
# 
# hrec %>% filter(kept==T) %>% group_by(population, chr, region) %>%
#   mutate(d_to_next_HR=lead(posl) - posr) %>%
#   group_by(population) %>%
#   summarise(d_to_next_HR=median(d_to_next_HR, na.rm=T))
# 
# 
# hrec %>% filter(kept==T) %>%
#   inner_join(carte %>% rename(poslw=posl, posrw=posr), by=c("chr","population","region")) %>%
#   filter(region !="C") %>%
#   group_by(population, chr, poslw, posrw) %>%
#   summarise(n=length(which(!(poslw >=posr | posrw <= posl)))) %>%
#   mutate(with_HR=ifelse(n>0,T, F)) %>%
#   ungroup() %>%
#   group_by(population) %>%
#   mutate(nb_HR=sum(n)) %>%
#   mutate(nb_w=n()) %>%
#   group_by(population, with_HR, nb_HR, nb_w) %>%
#   summarise(nb_HR_4Mb=quantile(n, 0.75), nb_w_HR=n()) %>%
#   filter(with_HR==T) %>%
#   mutate(prop_4Mb=nb_w_HR/nb_w, prop_tho=nb_HR/nb_w) 
# 
# 
# 
# #33CC99 = green
# #D7191C = blue
# #FF3333 = red
# #FF9900 = orange
# 
# 
# # 8) average historical rec rate 
# 
# 
# 
# fread(titre_resume) %>%
#   filter(w_center==T & population=="WE") %>%
#   na.omit() %>%
#   group_by(population) %>%
#   mutate(l=lpop*1e6) %>%
#   mutate(d=lambda_rho_med*l) %>%
#   mutate(region2=case_when(region %in% c("R1","R3")~"telo",
#                            region %in% c("R2a","R2b")~"peri",
#                            region=="C"~"centro")) %>%
#   group_by(region2) %>%
#   summarise(d=sum(d), l=sum(l), average_rho_WE=round((d/l)*1e3, digits=4)) %>%
#   dplyr::select(region2, average_rho_WE)
# 
# 
# fread(titre_resume) %>%
#   filter(w_center==T & population=="WE") %>%
#   na.omit() %>%
#   group_by(population) %>%
#   mutate(l=lpop*1e6) %>%
#   mutate(d=lambda_rho_med*l) %>%
#   mutate(genome=substr(chr,2,2)) %>%
#   group_by(genome) %>%
#   summarise(d=sum(d), l=sum(l), average_rho_WE=round((d/l)*1e3, digits=4)) %>%
#   dplyr::select(genome, average_rho_WE)
# 
# 
# fread(titre_resume) %>%
#   filter(w_center==T & region !="C") %>%
#   na.omit() %>%
#   group_by(population) %>%
#   mutate(l=lpop*1e6) %>%
#   mutate(d=lambda_rho_med*l) %>%
#   group_by(population) %>%
#   summarise(d=sum(d), l=sum(l), average_rho_pop=round((d/l)*1e3, digits=4)) %>%
#   dplyr::select(population, average_rho_pop)
# 
# 
# csre <- fread(titre_csre)
# 
# csre %>%
#   mutate(region2=case_when(region %in% c("R1","R3")~"telo",
#                            region %in% c("R2a","R2b")~"peri",
#                            region=="C"~"centro")) %>%
#   group_by(region2) %>%
#   summarise(l=sum(lpop), d=sum(dbay*1e2), average_rec_rate_csre=d/l) %>%
#   dplyr::select(region2,average_rec_rate_csre )
# 
# 
# 
# csre %>%
#   mutate(genome=substr(chr,2,2))%>%
#   group_by(genome) %>%
#   summarise(l=sum(lpop), d=sum(dbay*1e2), average_rec_rate_csre=d/l) %>%
#   dplyr::select(genome,average_rec_rate_csre )
# 
# 
# 
# rho <- fread(titre_resume) %>%
#   filter(w_center==T) %>%
#   na.omit() %>%
#   mutate(genome=substr(chr,2,2)) %>%
#   group_by(genome, population) %>%
#   mutate(l=lpop*1e6) %>%
#   mutate(d=lambda_rho_med*l)%>%
#   group_by(genome, population, region) %>%
#   summarise(d=sum(d), l=sum(l), average_rho=round((d/l)*1e3, digits=5)) %>%
#   dplyr::select(genome, population, region, average_rho) %>%
#   ungroup() %>%
#   mutate(region=factor(region, levels=c("R1","R3","R2a","R2b","C"))) %>%
#   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
#   spread(-population, value = "average_rho") %>%
#   ungroup()
# 
# 
# # tab2 <- fread(titre_csre) %>%
# #   mutate(genome=substr(chr,2,2)) %>%
# #   group_by(genome, region) %>%
# #   summarise(l=sum(lpop), d=sum(dbay*1e2), average_rec_rate_csre=round(d/l, digits=5)) %>%
# #   mutate(population="CsRe") %>%
# #   dplyr::select(genome, population, region, average_rec_rate_csre) %>%
# #   mutate(region=factor(region, levels=c("R1","R3","R2a","R2b","C"))) %>%
# #   spread(-population, value="average_rec_rate_csre") %>%
# #   ungroup() %>%
# #   rbind(rho,.) %>% arrange(genome, population)
# # 
# # write.table(tab2, titre_tab2, row.names = F, dec=".", sep="\t", quote=F)
# 
# 
# 
# 
# res <-   fread(titre_resume) %>%
#   filter(w_center==T) %>%
#   na.omit() 
# test <- fread(titre_domestication) %>%
#   rename(chr=V1, gene=V2, posl=V3,posr=V4, type=V5) %>%
#   mutate(taille=posr-posl) %>%
#   inner_join(res, by="chr") %>%
#   mutate(proxi=ifelse(!(posl >=posSNPrpop | posr <= posSNPlpop), T, F)) %>%
#   filter(proxi==T) %>%
#   dplyr::select(population, SNPlpop, SNPrpop, gene, type, lpop) %>%
#   unique() %>%
#   mutate(selection_signature=T) 
# 
# res2 <- res %>% full_join(test, by=c("population", "SNPlpop","SNPrpop", "lpop")) %>% 
#   dplyr::select(population, chr, posSNPlpop, posSNPrpop, lambda_rho_med, gene, type) %>%
#   inner_join(fread(titre_csre_4Mb) %>% dplyr::select(chr,posl,posr,cbay), by="chr") %>%
#   filter(!(posSNPlpop >= posr | posSNPrpop <= posl)) %>%
#   mutate(y=lambda_rho_med/(cbay*1e-8))
# 
# final <- data.frame()
# for (p in c("WE","EE","WA","EA")){
#   d <- res2 %>% filter(population==!!p)
#   mod <- rlm(log10(d$y)~1)
#   
#   proba_brutes <- pnorm(log10(d$y), mean=mod$coefficients, sd=mod$s, lower.tail = T)
#   proba_corrigees <- qvalue(proba_brutes)$qvalue
#   seuil_inf <- max(log10(d$y)[which(proba_corrigees<= 0.05)])
#   
#   
#   proba_brutes <- pnorm(log10(d$y), mean=mod$coefficients, sd=mod$s, lower.tail = F)
#   proba_corrigees <- qvalue(proba_brutes)$qvalue
#   seuil_sup <- min(log10(d$y)[which(proba_corrigees<= 0.05)])
#   
#   data <- data.frame(population=p, seuil_inf=seuil_inf, seuil_sup=seuil_sup)
#   
#   final <- rbind(final, data)
#   
# }
# 
# 
# graphe_dom <- res2 %>%
#   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
#   filter(is.na(type) | type=="domestication") %>%
#   ggplot(aes(x=log10(y))) +
#   geom_density(col="blue") +
#   geom_vline(data = final, aes(xintercept = seuil_inf)) +
#   geom_vline(data = final, aes(xintercept = seuil_sup)) +
#   geom_point(aes(x=log10(y), y=0, col=type)) +
#   scale_colour_manual("gene",values=c("red","blue"), na.translate=F) + 
#   geom_label_repel(aes(label = gene, y=0, col=type),
#                    box.padding   = 0.35, 
#                    point.padding = 0.5,
#                    segment.color = 'grey50')+
#   theme_light() +
#   xlab(expression(paste(rho,"/CsRe (log10)"))) +
#   facet_grid(population~.) +
#   theme(legend.position = "none") 
# #   
# # res %>%
# #   full_join(test, by=c("population", "SNPlpop","SNPrpop", "lpop")) %>%
# #   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
# #   filter(selection_signature==T & type=="domestication") %>%
# #   dplyr::select(population, lambda_med, gene) %>%
# #   spread(-gene, value = lambda_med) %>%
# #   dplyr::select(-gene) %>%
# #   cor(use="pairwise.complete.obs")
# 
# graphe_imp <- res2 %>%
#   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
#   filter(is.na(type) | type=="improvement") %>%
#   ggplot(aes(x=log10(y))) +
#   geom_density(col="blue") +
#   geom_vline(data = final, aes(xintercept = seuil_inf)) +
#   geom_vline(data = final, aes(xintercept = seuil_sup)) +
#   geom_point(aes(x=log10(y), y=0, col=type)) +
#   scale_colour_manual("gene",values=c("red","blue"), na.translate=F) + 
#   geom_label_repel(aes(label = gene, y=0, col=type),
#                    box.padding   = 0.35, 
#                    point.padding = 0.5,
#                    size=2,
#                    segment.color = 'grey50')+
#   theme_light() +
#   xlab(expression(paste(rho,"/CsRe (log10)"))) +
#   facet_grid(population~.) +
#   theme(legend.position = "none") 
# 
# graphe <- ggarrange(graphe_dom, graphe_imp, ncol=2, labels=c("A","B"))
# ggsave(titre_graphe_domestication, graphe)
# 
# 
# 
# 
# ##### table S1
# # tabs1 <- fread(titre_resume) %>%
# #   filter(region=="R3") %>%
# #   group_by(chr, region) %>%
# #   summarise(posr=round(max(posSNPrpop)/1e6)) %>%
# #   ungroup()
# # 
# # 
# # regions <- fread(titre_chr_tab, skip=1, header=T) 
# # 
# # tabs1 <- regions[,c(1:6)] %>%
# #   pivot_longer(cols=c("R1/R2a", "R2a/C","C/R2b","R2b/R3"), names_to = "region", values_to = "posr") %>%
# #   mutate(region=case_when(region=="R1/R2a"~"R1",
# #                           region=="R2a/C"~"R2a",
# #                           region=="C/R2b"~"C",
# #                           region=="R2b/R3"~"R2b")) %>%
# #   mutate(chr=substr(Chromosome, 4,5)) %>%
# #   dplyr::select(chr, region, posr) %>%
# #   rbind(tabs1,.) %>%
# #   mutate(region=factor(region, levels = c("R1","R2a","C","R2b","R3"))) %>%
# #   arrange(chr, region) %>%
# #   group_by(chr) %>%
# #   mutate(posl=lag(posr)) %>%
# #   mutate(posl=ifelse(is.na(posl),0,posl)) %>%
# #   mutate(size=posr-posl) %>%
# #   dplyr::select(chr, region, posl, posr, size)
# #   
# # tabs1 <- fread(titre_FST_blocs) %>%
# #   mutate(chr=substr(chrregion,1,2)) %>%
# #   mutate(region=substr(chrregion,3,5)) %>%
# #   dplyr::select(P1,P2, chr, region, FST) %>%
# #   full_join(tabs1, by=c("chr","region"))%>%
# #   dplyr::select(chr, region, P1, P2, posl, posr, size, FST) %>%
# #   rename(FST_haplotypes=FST)
# # 
# # 
# # 
# # tabs1 <- fread(titre_FST_SNP) %>%
# #   mutate(chr=substr(chrregion,1,2)) %>%
# #   mutate(region=substr(chrregion,3,5)) %>%
# #   dplyr::select(P1,P2, chr, region, FST) %>%
# #   full_join(tabs1, by=c("chr","region", "P1","P2"))%>%
# #   dplyr::select(chr, region, P1, P2, posl, posr, size, FST_haplotypes, FST) %>%
# #   rename(FST_SNP=FST) %>%
# #   mutate(P1=factor(P1, levels = c("WE","EE","WA",'EA'))) %>%
# #   mutate(P2=factor(P2, levels = c("WE","EE","WA",'EA')))%>%
# #   mutate(region=factor(region, levels = c("R1","R2a","C","R2b","R3"))) %>%
# #   arrange(chr, region, P1, P2) %>%
# #   mutate(FST_haplotypes=round(FST_haplotypes,2)) %>%
# #   mutate(FST_SNP=round(FST_SNP,2))
# # 
# # 
# # write.table(tabs1,titre_tab1, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
# # 
# 
# 
# ##### Annotation HR
# 
# 
# 
# 
# ##### Size of chromosomes
# 
# 
# g1 <- fread(titre_resume) %>%
#   filter(w_center==T) %>%
#   mutate(d=lambda_rho_med*(lpop*1e6)) %>%
#   group_by(chr, population) %>%
#   summarise(d=sum(d), l=sum(lpop), rec=d/(l*1e3)) %>%
#   mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
#   ggplot(aes(x=l, y=rec, col=population)) + geom_point(col="black") +
#   geom_smooth(method="lm", se=T) +
#   geom_text(aes(label=chr),hjust=0, vjust=0, col="black",size =5) +
#   theme_light() +
#   xlab("Physical size of chromosome (Mb)")+
#   ylab(expression(paste("Average ", rho, " (/kb)"))) +
#   facet_wrap(population~., ncol=4) +
#   scale_colour_manual(values=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA")) +
#   theme(legend.position = "top")
# 
# data <- fread(titre_resume) %>%
#   filter(w_center==T) %>%
#   mutate(d=lambda_rho_med*(lpop*1e6)) %>%
#   group_by(chr, population) %>%
#   summarise(d=sum(d), l=sum(lpop), rec=d/(l*1e3)) %>%
#   mutate(population=factor(population, levels=c("WE","EE","WA","EA")))
# mod <- lm(rec~l, data=data %>% filter(population=="WE"))
# summary(mod)
# mod <- lm(rec~l, data=data %>% filter(population=="EE"))
# summary(mod)
# mod <- lm(rec~l, data=data %>% filter(population=="WA"))
# summary(mod)
# mod <- lm(rec~l, data=data %>% filter(population=="EA"))
# summary(mod)
# 
# g2 <- fread(titre_csre) %>%
#   group_by(chr) %>%
#   summarise(d=sum(dbay*1e2), l=sum(lpop), rec=d/l) %>%
#   ggplot(aes(x=l, y=rec)) + geom_point()+
#   geom_smooth(method="lm", se=T) +
#   geom_text(aes(label=chr),hjust=0, vjust=0,size =5) +
#   theme_light() +
#   xlab("Physical size of chromosome (Mb)")+
#   ylab("Average CsRe Bayesian rec rates (cM/Mb")
# 
# data <- fread(titre_csre) %>%
#   group_by(chr) %>%
#   summarise(d=sum(dbay*1e2), l=sum(lpop), rec=d/l)
# mod <- lm(rec~l, data=data )
# summary(mod)
# 
# g <- ggarrange(g1,g2, nrow=2)
# ggsave(titre_graphe_chromosomes, g, height = 8, width = 20)
# 
# 
# 
# 
# ### lambda >=4 and lambda >=10 in two populations
# 
# 
# 
# cliques <- fread(titre_cliques_HR)
# # nb of cliques involving both WE and EA
# nWEEA <- cliques %>% filter(population%in% c("WE","EA") & coloc%in%c("WEEA", "WEEEEA","WEWAEA","WEEEWAEA")) %>% 
#   dplyr::select(clique) %>%
#   unique() 
# # nb cliques involving WE
# nWE <- cliques %>% filter(population%in% c("WE")) %>% 
#   dplyr::select(clique) %>%
#   unique() 
# # nb cliques involving EA
# nEA <- cliques %>% filter(population%in% c("EA")) %>% 
#   dplyr::select(clique) %>%
#   unique() 
# # nb cliques involving either WE or EA
# ntot <- rbind(nWE,nEA) %>% unique()
# length(which(nWEEA$clique %in% ntot$clique))/nrow(ntot)
# 
# nWEEE <- cliques %>% filter(population%in% c("WE","EE") & coloc%in%c("WEEE", "WEEEEA","WEEEWA","WEEEWAEA")) %>% 
#   dplyr::select(clique) %>%
#   unique() 
# # nb cliques involving WE
# nWE <- cliques %>% filter(population%in% c("WE")) %>% 
#   dplyr::select(clique) %>%
#   unique() 
# # nb cliques involving EA
# nEE <- cliques %>% filter(population%in% c("EE")) %>% 
#   dplyr::select(clique) %>%
#   unique() 
# # nb cliques involving either WE or EA
# ntot <- rbind(nWE,nEE) %>% unique()
# length(which(nWEEE$clique %in% ntot$clique))/nrow(ntot)
# 
# 
# final <- data.frame()
# 
# for (p in c("WE","EE","WA","EA")){
#   m <- cliques %>% filter(population==!!p) %>% summarise(m=max(lambda_med))
#   
#   for (i in seq(4,floor(m$m),1)) {
#     
#     partage <- cliques %>% filter(population==!!p & lambda_med >=i & coloc != !!p) %>% dplyr::select(IDintHR) %>%
#       unique()
#     total <- cliques %>% filter(population==!!p & lambda_med >=i) %>% dplyr::select(IDintHR) %>%
#       unique()
#     prop <- nrow(partage)/nrow(total)
#     
#     data=data.frame(population=p, seuil=i,prop=prop, partage=nrow(partage), total=nrow(total))
#     
#     final <- rbind(final, data)
#     
#   }
#   
# }
# 
# g <- ggplot(final, aes(x=seuil, y=prop, col=population)) + geom_line() + scale_x_continuous(trans="log10") +
#   xlab(expression(paste("Treshold for population's HR (", lambda,")"))) +
#   ylab(expression(paste("Proportion of HR (",lambda,"> treshold) overlapping HR(",lambda,">4) in other populations"))) +
#   theme_light() +
#   scale_colour_manual(values=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA"))
# 
# 
# ggsave(titre_graphe_HR_partages, g)
# 
# 
# final %>% filter(population=="WE") %>%
#   filter(seuil == 4 | seuil ==30)
# 
# 
sessionInfo()
# 
