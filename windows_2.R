
Sys.time()
cat("\n\nwindows_2.R\n\n")
rm(list = ls())
set.seed(1)
graphics.off()


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))







# Script qui définit une liste de fenetres, contenant des locus
# Fenetres chevauchantes aux extrémités seulement
# Fenetres de min 50 marqueurs, avec un bord supplémentaire d'au moins 1 marqueur
# ID des fenetres type : ID_cluster + num de la fenetre pour le cluster
# Fichier de synthèse en sortie



variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")









population <- variables[1]
chr <- variables[2]
dossier_map <- variables[3]
titre_pos <- variables[4]
dossier_sortie_f <- variables[5]
titre_fenetres_desc <- variables[6]


# population <-"EA"
# chr <- "1A"
# dossier_map <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/"
# titre_pos <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/SNP_positions.txt"
# dossier_sortie_f <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/fen/"
# titre_f <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/fen/fenetres.txt"


titre_map <- paste0(dossier_map,population,".map")
# titre_fenetres <- paste0(dossier_fenetres,population,"_",chr,".txt")

taille_min=5
taille_max_bord=30
taille_max_bord_cm=0.5
taille_max_fen=160

taille_max_centre=taille_max_fen - 2*taille_max_bord+2
taille_max_centre_cm=1
taille_min_fen=50
#taille_min_fen=100

pos <-  fread(titre_pos, header=T, dec=".", sep="\t", data.table = F)

m <- fread(titre_map, header=F, dec=".", sep="\t", data.table = F) %>%
  rename(chr_code_nombre=V1, SNP=V2, gen=V3, posSNP=V4) %>%
  left_join(pos, by="SNP", suffix=c("",".x")) %>%
  filter(chr %in% !!chr) %>%
  mutate(min=min(gen)) %>%
  mutate(gen=gen-min) %>%
  dplyr::select(chr, region, SNP, posSNP, gen)
  






limite <- function(m, nb_m, debut, taille_max, taille_cm){
  
  
  precedent <- m[debut, "gen"]
  
  
  diff_1 = taille_min
  
  
  diff_2 = min(which(m$gen > precedent + taille_cm), nb_m)- debut +1
  
  
  diff_3 <- taille_max
  
  
  
  if (diff_2 >= taille_min & diff_2 <= taille_max){
    
    
    out = min(which(m$gen > precedent + taille_cm), nb_m)
    
    
  } else if (diff_2 < taille_min){
    
    
    out = min(debut + taille_min -1, nb_m)
    
    
  } else if (diff_2 > taille_max){
    
    
    out = min(debut + taille_max -1, nb_m)
    
    
  }
  
  invisible(out)
  
  
}










limite_startw <- function(m, nb_m, startc){
  
  
  precedent <- m[startc, "gen"]
  
  diff_1 = taille_min
  
  
  diff_2 <- startc - max(which(m$gen < precedent - taille_max_bord_cm),1) +1
  
  
  diff_3 <- taille_max_bord
  
  
  
  if (diff_2 >= taille_min & diff_2 <= taille_max_bord){
    
    
    out = max(which(m$gen < precedent - taille_max_bord_cm),1)
    
    
  } else if (diff_2 < taille_min){
    
    
    out = max(startc - taille_min +1, 1)
    
    
  } else if (diff_2 > taille_max_bord){
    
    
    out = min(startc - taille_max_bord +1, nb_m)
    
    
  }
  
  invisible(out)
  
  
}



fenetre <- function(m, nb_m, startc,f){
  
  startw <- limite_startw(m, nb_m, startc)
  
  
  stopc <- limite(m, nb_m, startc, taille_max_centre, taille_max_centre_cm)
  
  stopw <- limite(m, nb_m, stopc, taille_max_bord, taille_max_bord_cm)
  
  
  
  
  
  # Cas où stopc == stopw
  
  if(stopc == stopw & stopc - startc +1 >= 4){
    
    stopc = stopc -2
    
    
  } 
  
  
  
  
  
  
  while(stopw - startw +1 < taille_min_fen & stopw < nb_m){
    
    k=0
    
    
    if (startc - startw +1 < taille_max_bord & startw>1){
      
      startw <- max(startw -1, 1)
      
      k=1
    }
    
    
    if (stopw - startw +1 < taille_min_fen & stopw - stopc +1 < taille_max_bord & k ==0 & stopw < nb_m){
      
      
      stopw <- min(stopw +1, nb_m)
      
      
      
      k=1
    }
    
    
    if (k==0) {
      
      
      if (sample(c("l","r"), 1)=="l"){
        
        startw <- max(startw -1, 1)
        
      } else {
        stopw <- min(stopw +1, nb_m)
      }
      
    }
    
    
    
    
    
    
    
  }
  
  if (f==1 & startw >1){
    
    startw2 <- 1
    
    startc2 <- startc - startw +1
    
    
    stopc2 <- stopc - startw +1
    
    stopw2 <- stopw - startw +1
    
    startw <- startw2
    stopw <- stopw2
    startc <- startc2
    stopc <- stopc2
    
  }
  
  
  
  poslfen=m$pos[startw]
  posrfen=m$pos[stopw]
  centerlfen=m$pos[startc]
  centerrfen=m$pos[stopc]
  
  genl=m$gen[startw]
  genr=m$gen[stopw]
  center_genl=m$gen[startc]
  center_genr=m$gen[stopc]
  
  nb_mqs = stopw - startw +1
  
  ID=paste0(population,"_", chr,"_",f)
  
  region=m$region[startw]
  
  
  fen <- data.frame(population,
                    chr,
                    poslfen,
                    posrfen,
                    centerlfen,
                    centerrfen,
                    genl,
                    genr,
                    center_genl,
                    center_genr,
                    nb_mqs,
                    startw,
                    startc,
                    stopc,
                    stopw,
                    ID,
                    region)
  
  invisible(fen)
  
  
  
}






nb_m <- nrow(m)
startc <- min(which(m$gen > taille_max_bord_cm))
stopw <- 0
f=0
fenf <- data.frame()
while(stopw < nb_m - 2*taille_min){
  #while(f < 162-1 ){
  
  f=f+1
  
  
  fen <- fenetre(m, nb_m, startc, f)
  
  fenf <- rbind(fenf, fen)
  
  startc <- max(fenf$stopc)
  stopw <- max(fenf$stopw)
  
  
  print(f)
  
  titre_f=paste0(dossier_sortie_f,"w_",fen$ID,".txt")
  
  
  mqs <- m[which(m$posSNP >= fen$poslfen & m$posSNP <= fen$posrfen),"SNP"]
  
  write.table(mqs, titre_f, col.names = F, row.names = F, dec=".", sep="\t", quote=F)
  
  
}


# Pour prendre tous les mqs du chromosome
 while (last(fenf$nb_mqs) < taille_max_fen & max(fenf$stopw) < nb_m){
   
   fenf[f, "stopw"] <- max(fenf$stopw) +1
   
 }
 
 stopw <- fenf[f, "stopw"]
 fenf[f,"posrfen"] <- m[stopw,"posSNP"]
 fenf[f,"genr"] <- m[stopw,"gen"]
 fenf[f,"nb_mqs"] <- stopw - fen$startw +1


mqs <- m[which(m$posSNP >= max(fen$poslfen) & m$posSNP <= max(fen$posrfen)),"SNP"]

write.table(mqs, titre_f, col.names = F, row.names = F, dec=".", sep="\t", quote=F)


fenf$taillel <- fenf$startc - fenf$startw +1
fenf$tailler <- fenf$stopw - fenf$stopc +1
fenf$taillec <- fenf$stopc - fenf$startc +1
fenf$taille <- fenf$stopw - fenf$startw +1
fenf$taillegencenter <- fenf$center_genr - fenf$center_genl
fenf$taillegenl <- fenf$center_genl - fenf$genl
fenf$taillegenr <- fenf$genr - fenf$center_genr
fenf$taillephyl <- fenf$centerlfen - fenf$poslfen
fenf$taillephyr <- fenf$posrfen - fenf$centerrfen
fenf$taillephycenter <- fenf$centerrfen - fenf$centerlfen
fenf$taillephytotal <- fenf$posrfen - fenf$poslfen


cat("\n\n nb mqs left : \n")
quantile(fenf$taillel)
cat("\n\n nb mqs right : \n ")
quantile(fenf$tailler)
cat("\n\n nb mqs center : \n")
quantile(fenf$taillec)
cat("\n\n nb mqs total : \n")
quantile(fenf$taille)
cat("\n\n genetique size left : \n")
quantile(fenf$taillegenl)
cat("\n\n genetique size right : \n")
quantile(fenf$taillegenr)
cat("\n\n genetique size center : \n")
quantile(fenf$taillegencenter)
cat("\n\n")
cat("\n\n physical size left : \n")
quantile(fenf$taillephyl)
cat("\n\n")
cat("\n\n physical size right : \n")
quantile(fenf$taillephyr)
cat("\n\n")
cat("\n\n physical size center : \n")
quantile(fenf$taillephycenter)
cat("\n\n")
cat("\n\n total physical size : \n")
quantile(fenf$taillephytotal)
cat("\n\n")



fenf <- fenf %>%
  dplyr::select(population,chr,region, poslfen,posrfen,centerlfen,centerrfen,nb_mqs,ID)


head(fenf)
tail(fenf)


write.table(fenf, titre_fenetres_desc, col.names = T, row.names = F, dec=".", sep="\t", quote=F)


sessionInfo()





# head(fenf)
# max(fenf$taille)
# max(fenf$taillec)
# max(fenf$tailler)
# max(fenf$taillel)
# min(fenf$taillel)
# min(fenf$tailler)
# min(fenf$taillec)
# min(fenf$taille)

