
Sys.time()
rm(list = ls())
graphics.off()
set.seed(1)
variables <- commandArgs(trailingOnly=TRUE)


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))




cat("\n\ninterpolation.R\n\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")



titre_csre_bayesian <-variables[1]
titre_pos_landraces <-variables[2]
titre_dossier_png <- variables[3]
titre_update_pos <-  variables[4]
titre_update_gen <-  variables[5]
titre_SNP_PLINK <-  variables[6]
titre_graphe <-  variables[7]




 # titre_csre_bayesian <-"/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_genetic_map.txt"
 # titre_pos_landraces <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/SNP_positions.txt"
 # titre_dossier_png <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/interpolation/"
 # titre_update_pos <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/update_pos.txt"
 # titre_update_gen <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/update_gen.txt"
 # titre_SNP_PLINK <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/SNP_PLINK.txt"
 # titre_graphe <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/graphes/bayesian_and_interpolation.png"

#  titre_csre_bayesian <- "/work/adanguy/these/pipeline/160320/sources/csre_classique.txt"       
#  titre_pos <- "/work/adanguy/these/pipeline/160320/sources/SNP_positions.txt" 
#  titre_sortie <- "/work/adanguy/these/pipeline/030420/sources/SNP_positions.txt" 
#  titre_dossier_png <- "/work/adanguy/these/pipeline/030420/graphes/interpolation/"    
#  titre_update_pos <- "/work/adanguy/these/pipeline/030420/temp/autres/update_pos.txt"
#  titre_update_gen <- "/work/adanguy/these/pipeline/030420/temp/autres/update_gen.txt"
#  titre_SNP_PLINK <- "/work/adanguy/these/pipeline/030420/temp/autres/SNP_PLINK.txt" 




# dataframe with mapped markers
cat("\n\n Input 1 : Empirical and bayesian estimates of meiotic recombinnation rates of CsRe \n\n")
csre <- fread(titre_csre_bayesian, header=T, dec=".", sep="\t", data.table = F) 
head(csre) %>% as.data.frame()
csre <- csre %>%
  arrange(chr, posSNPlpop) %>%
  mutate(d=cbay*lpop) %>%   # calcul de la distance genetique pour chaque intervalle
  group_by(chr) %>%
  mutate(d_cumul=cumsum(d)) %>% # calcul de la distance genetique cumulee a chaque intervall
  mutate(posSNP=posSNPrpop) %>%
  dplyr::select(chr, region, posSNP, d_cumul)

csre2 <-  fread(titre_csre_bayesian, header=T, dec=".", sep="\t", data.table = F) %>%
  arrange(chr, posSNPlpop) %>%
  group_by(chr) %>%
  slice(1) %>%
  mutate(posSNP=posSNPlpop, d_cumul=0) %>%
  dplyr::select(chr, region, posSNP, d_cumul) %>%
  rbind(., csre) %>%
  arrange(chr, posSNP)

# dataframe of data to interpolate
cat("\n\n Input 2 : physical position of SNP for interpolation \n\n")
pos <- fread(titre_pos_landraces, header=T, dec=".", sep="\t", data.table = F) 
head(pos)
pos <- pos %>%
  dplyr::select(chr, chr_code_nombre, region, SNP, posSNP, other_SNP_ID, category) %>%
  full_join(csre2, by=c("chr", "posSNP", "region")) %>% # possibilite d'avoir des mqs differents dans les 2 jeux de donnees
  dplyr::select(chr, region, posSNP, d_cumul) %>%
  arrange(chr, posSNP) %>%
  ungroup()





na_index <- pos %>%
  group_by(chr) %>%
  arrange(posSNP) %>%
  mutate(interpolation=approxfun(x=posSNP, y=d_cumul, ties="ordered", rule=1)(posSNP)) %>% # Fonction interpolation lineaire
  summarise(first_non_NA=first(which(!is.na(interpolation))), last_non_NA=last(which(!is.na(interpolation)))) %>%
  ungroup()

# na_index ne sert qu'a obtenir les mqs pour lesquels on peut pas calculer d'interpolation
# Ces mqs sont positionnes posSNPsiquement avant les mqs de csre (en R1) ou après (en R3) où interpolation linéaire n'est possible
# first_non_NA = premier mqs a avoir une d_cumul != NA
# last_non_NA = dernier mqs a avoir une d_cumul != NA
# first et last sont calcules pour chaque chromosome
# Les mqs au dela de first_non_NA et en dessous de last_non_NA auront une position génétique interpolee linéairement
# les autres mqs auront une position génétique qui sera proportionelle à leur distance posSNPsique et au taux de rec moyen de R1 ou R3

m <- data.frame()
chr="1A"
# On va traiter différement les portions des chromosomes left et right
# Les deux portions auront la mm interpolations pour les mqs situés entre ceux de csre
# mais les modeles lineaires pour estimer les positions des mqs situés avant ou après csre changent
for (chr in unique(pos$chr)){
  
  print(chr)
  
  
  if(na_index$first_non_NA[which(na_index$chr==chr)] ==1){ # si tous les mqs en R1 ont pas pu etre interpoles (first_non_NA=1)
    
    # Dans ce cas, il suffit d'utiliser l'interpolation
    
    left <- pos %>%
      filter(chr==!!chr) %>%
      mutate(interpolation=approxfun(x=posSNP, y=d_cumul, ties="ordered", rule=1)(posSNP)) %>%
      mutate(region = factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
      dplyr::select(chr, region, posSNP, d_cumul, interpolation)
    
    
  } else { # si tous les mqs en R1 ont pas pu etre interpoles (first_non_NA !=1)
    
    # Dans ce cas, il faut estimer le taux de rec moyen de R1
    
    
    mod_left <- pos %>% 
      filter(chr==!!chr) %>%
      filter(region=="R1") %>%
      filter(!is.na(d_cumul)) %>%
      lm(data=., d_cumul~posSNP) # taux de rec moyen de R1
    
    mod_left$coefficients[1] <- 0
    
    left <- pos %>%
      filter(chr==!!chr) %>%
      mutate(interpolation=approxfun(x=posSNP, y=d_cumul, ties="ordered", rule=1)(posSNP)) %>%
      left_join(na_index, by="chr") %>%
      arrange(posSNP) %>%
      mutate(interpolation2=ifelse(row_number() <= first_non_NA,
                                   predict(mod_left, newdata=data.frame(posSNP)), # modele lineaire pour les mqs extreme. On enleve l'intercept pour avoir juste le taux moyen
                                   interpolation)) %>% # interpolation pour les mqs au centre du chr
      mutate(jointure=interpolation2 - interpolation) %>% # l'interpolation et le modele lineaire ne prédisent pas necessairement la mm valeur au niveau de la joinnture left-centre du chr.
      mutate(jointure=ifelse(jointure ==0, NA, jointure)) %>%
      group_by(chr) %>%
      mutate(jointure=max(jointure, na.rm=T)) %>%
      ungroup() %>%
      mutate(interpolation3=case_when(row_number() > first_non_NA~interpolation2,
                                      row_number() <=first_non_NA & jointure >=0~interpolation2-jointure,
                                      row_number() <=first_non_NA & jointure <0~interpolation2 + jointure)) %>%
      mutate(minimum=interpolation3[1]) %>%
      mutate(interpolation4=ifelse(minimum <= 0, interpolation3 + abs(minimum), interpolation3 - minimum )) %>% # on fait commencer la carte genetique à 0 cM
      mutate(interpolation5=ifelse(row_number() > last_non_NA, NA, interpolation4)) %>%
      mutate(interpolation=interpolation5) %>%
      ungroup() %>%
      mutate(region = factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
      dplyr::select(chr, region, posSNP, d_cumul, interpolation)
    
    
    
  }
  
  
  if(na_index$last_non_NA[which(na_index$chr==chr)] ==nrow(left)){ # si tous les mqs en R1 ont pas pu etre interpoles (last_non_NA=dernier mqs)
    
    # Pas besoin d'en faire plus
    
    right <- left
    
  } else { # si tous les mqs en R3 ont pas pu etre interpoles 
    
    # Dans ce cas, il faut estimer le taux de rec moyen de R3
    
    
    
    mod_right <- pos %>%
      filter(chr==!!chr) %>%
      filter(region=="R3") %>%
      filter(!is.na(d_cumul)) %>%
      lm(data=., d_cumul~posSNP)# taux de rec moyen de R3
    
    
    mod_right$coefficients[1] <- 0
    
    
    
    right <-  left %>% left_join(na_index, by="chr") %>%
      arrange(posSNP) %>%
      mutate(interpolation2=ifelse(row_number() >= last_non_NA,
                                   predict(mod_right, newdata=data.frame(posSNP)) ,
                                   interpolation)) %>%
      mutate(jointure=interpolation2 - interpolation) %>%
      group_by(chr) %>%
      mutate(jointure=ifelse(jointure !=0, jointure, NA)) %>%
      mutate(jointure=max(jointure, na.rm=T)) %>%
      ungroup() %>%
      mutate(interpolation3=case_when(row_number() < last_non_NA~interpolation2,
                                      row_number() >=last_non_NA & jointure >=0~interpolation2 - jointure,
                                      row_number() >=last_non_NA & jointure <0~interpolation2 - jointure))%>%
      mutate(interpolation=interpolation3) %>%
      ungroup() %>%
      mutate(region = factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
      dplyr::select(chr, region, posSNP, d_cumul, interpolation)
    
    
  } 
  
  
  
  
  
  
  
  
  
  titre_png <- paste0(titre_dossier_png,"interpolation", chr,".png")
  png(titre_png)
  
  nb_mqs <- nrow(right)
  nb_mqs_gen <- right %>% filter(!is.na(d_cumul)) %>% ungroup() %>% count()
  position_legend <- max(right$posSNP)
  
  
  graphe <- right %>% ggplot(aes(x=posSNP, y=interpolation)) +
    geom_line(col="black") +
    geom_point(aes(x=posSNP, y=d_cumul, col=region)) +
    theme_light() +
    xlab("physical position (pb)") +
    ylab("genetic position (cM)") +
    ggtitle(chr) +
    guides(colour=guide_legend("CsRe SNPs")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x=position_legend, y = 10, vjust=1, hjust=1, label =  c(paste0("CsRe ", nb_mqs_gen," SNPs \n Interpolation ", nb_mqs," SNPs")))
  
  if (chr=="1A"){
  cat("\n\n Graph 1 : Graphical representation of interpolation \n\n")
  }
  print(graphe)
  dev.off()
  
  m <- rbind(m, right)
  
  
} 


head(m)



pos_landraces <- fread(titre_pos_landraces, header=T, dec=".", sep="\t", data.table = F) %>%
  dplyr::select(chr,chr_code_nombre, region, SNP, posSNP,other_SNP_ID,category) %>%
  left_join(m, by=c("chr", "posSNP","region"))%>%
  dplyr::select(chr,chr_code_nombre, region, SNP, posSNP,other_SNP_ID,category, d_cumul, interpolation) 





cat("\n\n Output 1 : Interpolated genetic position for SNP \n\n")
write.table(pos_landraces, titre_pos_landraces, sep="\t", col.names = T, row.names = F, dec=".", quote=F)
head(pos_landraces)


update_pos <- pos_landraces %>% dplyr::select(SNP, posSNP)
update_gen <- pos_landraces %>% dplyr::select(SNP, interpolation)
SNP_PLINK <- pos_landraces %>% dplyr::select(SNP)

cat("\n\n Output 3 : updated physical position for the PLINK files \n\n")
write.table(update_pos, titre_update_pos, sep="\t", col.names = F, row.names = F, dec=".", quote=F)
head(update_pos)
cat("\n\n Output 4 : updated genetic position for the PLINK files \n\n")
write.table(update_gen, titre_update_gen, sep="\t", col.names = F, row.names = F, dec=".", quote=F)
head(update_gen)
cat("\n\n Output 5 : list of SNP to keep in PLINK files \n\n")
write.table(SNP_PLINK, titre_SNP_PLINK, sep="\t", col.names = F, row.names = F, dec=".", quote=F)
head(SNP_PLINK)

# 
# i=0
# tab <- data.frame()
# size=1*1e6
# step=1*1e6
# t <-  pos_landraces  %>% filter(chr=="3B")
# sup <- 0
# while(sup <= max(t$posSNP)) {
#   i=i+1
#   
#   inf <-  step*(i-1)
#   sup <- step*(i-1) + size
#   
#   temp <- t %>% filter(posSNP >= inf & posSNP <= sup)
#   
#   if (nrow(temp) >0){
#     dmin <- min(temp$interpolation)
#     dmax <- max(temp$interpolation)
#     lmin <- min(temp$posSNP)
#     lmax <- max(temp$posSNP)
#     region <- as.data.frame(table(temp$region)) %>% arrange(desc(Freq)) %>% slice(1) %>% dplyr::select(Var1) %>% unlist() %>% as.vector()
#     
#     rec <- (dmax-dmin)/((lmax-lmin)/1e6)
#     
#     tab[i,"sup"] <- sup
#     tab[i, "rec"] <- rec
#     tab[i,"region"] <- region
#   }
# }
# 

sessionInfo()

