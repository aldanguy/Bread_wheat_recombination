
# module load system/R-3.4.3_bis


### Compute FST based on Weir & Cockerham (1984) distance of genomic regions based on haplotypic alleles



Sys.time()
cat("\n\nFST_blocks.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)

variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(hierfstat))
suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(corrplot))





titre_regions <- variables[1]
titre_correspondance_chr <-variables[2]
titre_haplo <-variables[3]
titre_bornes <- variables[4]
titre_landraces <-variables[5]
titre_stat_chrregion <- variables[6]
titre_matrice_Fst <- variables[7]



# titre_regions <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Decoupage_chr_ble.tab"
# titre_correspondance_chr <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Codes_chr.txt"
#    titre_haplo <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/landrace632_8741hap_sophie.var"
#    titre_landraces<- "/home/adanguydesd/Documents/These_Alice/pipeline/sources/landraces.txt"
#    titre_bornes <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/stat_haplo8741.txt"
#    titre_sortie <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/FST.txt"
# 
         # titre_regions <- "/work/adanguy/these/pipeline/amont/Decoupage_chr_ble.tab"            
         # titre_correspondance_chr <- "/work/adanguy/these/pipeline/amont/Codes_chr.txt"                    
         # titre_haplo <- "/work/adanguy/these/pipeline/amont/landrace632_8741hap_sophie.var"   
         # titre_bornes <- "/work/adanguy/these/pipeline/amont/stat_haplo8741.txt"               
         # titre_landraces <- "/work/adanguy/these/pipeline/030420/sources/landraces.txt"           
         # titre_stat_chrregion <- "/work/adanguy/these/pipeline/030420/FST/chrregion_FST_haplotypes.txt"
         # titre_stat_bloc <- "/work/adanguy/these/pipeline/030420/FST/blocs_FST_haplotypes.txt"    
         # titre_stat_pop <- "/work/adanguy/these/pipeline/030420/FST/blocs_he_haplotypes.txt"  
         # titre_matrice_Fst <- "/work/adanguy/these/pipeline/030420/FST/matrice_Fst_blocs.txt"  

correspondance_chr <- read.table(titre_correspondance_chr, header=F, dec=".", sep="\t") %>%
  rename( chr = V1, chr_code_nombre = V2)




regions  <- read.table(titre_regions, header=T, dec=".", sep="\t", skip=1) %>%
  dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
  mutate(chr = str_remove(Chromosome,"chr")) %>%
  dplyr::select(-one_of("Chromosome")) %>% 
  full_join(correspondance_chr, by="chr") %>%
  rename(R1=R1.R2a) %>%
  rename(R2a=R2a.C) %>%
  rename(C=C.R2b) %>%
  rename(R2b=R2b.R3) %>%
  mutate(R3=1e9) %>%
  mutate(R1=R1*1e6) %>%
  mutate(R2a=R2a*1e6) %>%
  mutate(C=C*1e6) %>%
  mutate(R2b=R2b*1e6) %>%
  pivot_longer(-c("chr_code_nombre","chr"), names_to ="region", values_to = "pos_max") %>%
  group_by(chr_code_nombre) %>%
  mutate(pos_min=lag(pos_max) +1) %>%
  mutate(pos_min=ifelse(!is.na(pos_min),pos_min,1)) %>%
  mutate(chrregion=paste0(chr,region)) %>%
  mutate(pos=pos_min+((pos_max-pos_min)/2)) %>%
  dplyr::select(chr, chr_code_nombre, chrregion, region, pos_min, pos_max) 





b <- fread(titre_bornes, header=T, dec=".", sep="\t") %>%
  rename(chr_code_nombre=chr) %>%
  full_join(correspondance_chr, by=c("chr_code_nombre")) %>%
  mutate(bloc=paste0("B", seq(1,nrow(.),1)))

haplo <- fread(titre_haplo, header=F, dec=".", sep="\t", colClasses = "character")
colnames(haplo) <- c("Unit", paste0(rep(paste0("B", 1:((ncol(haplo)-1)/2)), each=2),".",1:2))



landraces <- fread(titre_landraces, header=T, dec=".", sep="\t") %>%
  arrange(Unit) 

l <- fread(titre_landraces, header=T, dec=".", sep="\t") %>%
  arrange(Unit)%>%
  dplyr::select(Unit, ERGE) %>%
  mutate(Unit = as.character(Unit)) %>%
  full_join(haplo, by="Unit") %>%
  dplyr::select(-Unit) 

# l <- l[,c(1,1:361)]

l <- l %>% mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("999"), replacement = "NA"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("1"), replacement = "01"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("2"), replacement = "02"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("3"), replacement = "03"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("4"), replacement = "04"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("5"), replacement = "05"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("6"), replacement = "06"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("7"), replacement = "07"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("8"), replacement = "08"))) %>%
  mutate_at(vars(starts_with("B")), funs(str_replace_all(., pattern = c("9"), replacement = "09"))) %>%
  mutate_at(vars(starts_with("B")), funs(as.numeric(.)))
  





ERGE <- landraces %>% filter(kept==T) %>%
  dplyr::select(ERGE) %>% 
  unlist() %>%
  as.vector()


# ERGE <- landraces %>%
#   dplyr::select(ERGE) %>% 
#   unlist() %>%
#   as.vector()


population <- landraces %>% filter(kept==T) %>%
  dplyr::select(population) %>% 
  unlist() %>%
  as.vector()

# population <- landraces  %>%
#   dplyr::select(STRgp4) %>% 
#   unlist() %>%
#   as.vector()

l <- l[which(l$ERGE %in% ERGE),]
l$population <- population
l <- l[,c(ncol(l),c(1:(ncol(l)-1)))]
l <- l[,c(1,4:ncol(l))]
l[1:5,1:11]


# 
# # Global paiwise FST
# tab
# 
matrice1all <- pairwise.WCfst(l, diploid=T)
diag(matrice1all) <- 0
cat("\n\n Weir & Cockerham\n")
genet.dist(l,diploid=T, method="WC84")
cat("\n\n Reynolds\n")
genet.dist(l,diploid=T, method="Fst")
cat("\n\n Nei\n")
genet.dist(l,diploid=T, method="Nei87")


# 

col.ordre <- c("WE","EE","WA","EA")
matrice1all <- matrice1all[col.ordre,col.ordre]
#matrice1all

# png(titre_graphe, width = 6, height = 6, units = 'in', res = 300)
# 
# col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
#                            "#33FFCC", "#007FFF", "#6699FF"))
# corrplot(matrice1all, method="color", col=col1(10),  
#          type="upper", 
#          addCoef.col = "black", # Ajout du coefficient de corrélation
#          tl.col="black", tl.srt=45, sig.level = 0.01, insig = "blank", 
#          # Cacher les coefficients de corrélation sur la diagonale
#          diag=FALSE , is.corr = F, shade.lwd=0,number.digits=3
# )
# 
# dev.off()

stat_chrregion <- data.frame()


c="2AC"

# Calcul du FST par paire de pop dans les chr-region
for (c in unique(regions$chrregion)){
  
  print(c)
  
  pmin <- regions$pos_min[which(regions$chrregion==c)]
  pmax <- regions$pos_max[which(regions$chrregion==c)]
  
  chr=substr(c, 1,2)
  
  
  blocs <- b %>% filter(chr==!!chr) %>% filter(!(!!pmin >= stop | !!pmax <= start)) %>% dplyr::select(bloc) %>% unlist() %>% as.vector()
  # blocs <- b %>% filter(chr==!!chr) %>% filter(start <= !!pmax) %>% dplyr::select(bloc) %>% unlist() %>% as.vector()
  
  blocs2 <- colnames(l)[which(colnames(l) %in% paste0(rep(blocs, each=2), ".", 1:2))]
  

  temp2 <- l %>% dplyr::select(population, one_of(blocs2))
  

  
  nb_blocs=length(blocs)
  print(nb_blocs)
  
  if (nb_blocs >=1){
    
    
    
    
    #FST
    matrice1 <- pairwise.WCfst(temp2, diploid=T)
    diag(matrice1) <- 0
    matrice2 <- matrice1[lower.tri(matrice1)]
    stat_chrregion_temp <- data.frame(P1=t(combn(row.names(matrice1), 2))[,1],
                                      P2=t(combn(row.names(matrice1), 2))[,2],
                                      chrregion=rep(c, times=nrow(t(combn(row.names(matrice1), 2)))),
                                      pos=pmin,
                                      nb_blocs=nb_blocs,
                                      FST=matrice2)
    
    
    
  } else {
    
    stat_chrregion_temp <- data.frame(P1=t(combn(row.names(matrice1), 2))[,1],
                                      P2=t(combn(row.names(matrice1), 2))[,2],
                                      chrregion=rep(c, times=nrow(t(combn(row.names(matrice1), 2)))),
                                      pos=pmin,
                                      nb_blocs=0,
                                      FST=rep(NA, times=nrow(t(combn(row.names(matrice1), 2)))))
    
    
  }
  
  stat_chrregion <- rbind(stat_chrregion, stat_chrregion_temp)
  
  
}





stat_chrregion <- stat_chrregion %>%
  mutate(P1.2 = case_when(P1=="WE"|P2=="WE" ~ "WE",
                          (P1!= "WE" & P2 !="WE") & (P1=="EE" | P2 =="EE") ~ "EE",
                          (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE") & (P1=="WA"|P2=="WA") ~ "WA",
                          (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA") & (P2=="EA"|P2=="EA")~ "EA")) %>%
  mutate(P2.2 = ifelse(P1==P1.2, as.character(P2), as.character(P1))) %>%
  dplyr::select(-one_of("P1","P2")) %>%
  rename(P1=P1.2, P2=P2.2) %>%
  mutate(methode="haplotype") %>%
  dplyr::select(chrregion, P1, P2, FST, nb_blocs, methode)


cat("\n\n FST for each chrregion \n")
head(stat_chrregion)
write.table(stat_chrregion, titre_stat_chrregion, col.names=T, row.names=F, dec=".", sep="\t", quote=F)

cat("\n\n FST matrix \n")
matrice1all
write.table(matrice1all, titre_matrice_Fst, col.names=T, row.names=T, dec=".", sep="\t", quote=F)





# 
# 
# #### Per haplotypic bloc
# stat_bloc <- data.frame()
# stat_pop <- data.frame()
# 
# 
# 
# # Calcul du FST par paire de pop dans les blocs et calcul du he par bloc et par pop
# for (i in 1:nrow(b)){
#   
#   print(i)
#   
#   
#   blocs1 <- b$bloc[i]
#   blocs2 <- paste0(rep(b$bloc[i], each=2), ".", 1:2)
#   
#   if (length(which(colnames(l2) %in% blocs1)) >0){ # utile juste pour quand je fais des tests
#     
#     temp <- l2 %>% dplyr::select(one_of(blocs1))
#     
#     temp2 <- l %>% dplyr::select(one_of(blocs2))
#     
#     
#     
#     alleles=unique(c(as.character(temp2[,1]), as.character(temp2[,2])))
#     nb_alleles = length(alleles[-min(which(alleles == "NA"),1e9)])
#     
#     
#     
#     stat_pop_temp <- data.frame()
#     k=0
#     
#     for (p in unique(population)){
#       k=k+1
#       
#       
#       tempp <- as.matrix(temp[which(population == p),])
#       
#       temp2p <- temp2[which(population == p),]
#       
#       
#       alleles=unique(c(as.character(temp2p[,1]), as.character(temp2p[,2])))
#       nb_alleles = length(alleles[-min(which(alleles == "NA"),1e9)])
#       tab <- df2genind(tempp, ploidy=2, sep=" ", ncode=4, NA.char = "NA")
#       
#       if (nb_alleles >1){
#         
#         he <- Hs(tab)
#         
#       } else {
#         
#         he <- NA
#       }
#       
#       stat_pop_temp[k,"population"] <- p
#       stat_pop_temp[k,"bloc"] <- b$bloc[i]
#       stat_pop_temp[k,"he"] <- he
#       stat_pop_temp[k,"pos"] <- b$pos[i]
#       stat_pop_temp[k,"nb_alleles"] <- nb_alleles
#       stat_pop_temp[k,"methode"] <- "haplotype"
#       
#     }
#     
#     
#     if(length(which(complete.cases(stat_pop_temp)))>=3){
#       
#       #FST
#       tab <- df2genind(temp, ploidy=2, sep=" ", ncode=4, pop=population, NA.char = "NA")
#       
#       matrice1 <- pairwise.fst(tab, res="matrix")
#       matrice2 <- matrice1[lower.tri(matrice1)]
#       stat_bloc_temp <- data.frame(P1=t(combn(row.names(matrice1), 2))[,1],
#                                    P2=t(combn(row.names(matrice1), 2))[,2],
#                                    bloc=b$bloc[i],
#                                    pos=b$pos[i],
#                                    nb_alleles=nb_alleles,
#                                    FST=matrice2,
#                                    mutate(methode="haplotype")) 
#       
#       
#     }else if(length(which(complete.cases(stat_pop_temp)))<=2 & length(which(complete.cases(stat_pop_temp))) >=1){  # Si 1 à 3 populations sont monomorphic, le pairwise.FST ne fonctionne plus
#       
#       pop_a_pb <- stat_pop_temp %>% filter(is.na(he)) %>% dplyr::select(population) %>% unlist() %>% as.vector()
#       pop_pas_a_pb <- stat_pop_temp %>% filter(!is.na(he)) %>% dplyr::select(population) %>% unlist() %>% as.vector()
#       
#       
#       stat_bloc_temp1 <- data.frame()
#       
#       for (c in pop_a_pb){
#         
#         tempp <- as.matrix(temp[which(population %in% c(c,pop_pas_a_pb)),])
#         
#         temp2p <- temp2[which( population %in%c(c,pop_pas_a_pb)),]
#         
#         tab <- df2genind(tempp, ploidy=2, sep=" ", ncode=4, pop=population[which( population %in%c(c,pop_pas_a_pb))], NA.char = "NA")
#         
#         matrice1 <- pairwise.fst(tab, res="matrix")
#         matrice2 <- matrice1[lower.tri(matrice1)]
#         stat_bloc_temp2 <- data.frame(P1=t(combn(row.names(matrice1), 2))[,1],
#                                       P2=t(combn(row.names(matrice1), 2))[,2],
#                                       bloc=b$bloc[i],
#                                       pos=b$pos[i],
#                                       nb_alleles=nb_alleles,
#                                       FST=matrice2,
#                                       mutate(methode="haplotype"))
#         
#         
#         stat_bloc_temp1 <- rbind(stat_bloc_temp1,stat_bloc_temp2) 
#         
#       }
#       
#       
#       
#       stat_bloc_temp3 <- data.frame(P1=t(combn(pop_a_pb, 2))[,1],
#                                     P2=t(combn(pop_a_pb, 2))[,2],
#                                     bloc=b$bloc[i],
#                                     pos=b$pos[i],
#                                     nb_alleles=nb_alleles,
#                                     FST=NA,
#                                     mutate(methode="haplotype")) 
#       
#       stat_bloc_temp <- rbind(stat_bloc_temp1,stat_bloc_temp3)  %>%
#         unique()  
#       
#       
#       
#       
#       
#       
#     }
#     
#     stat_bloc <- rbind(stat_bloc, stat_bloc_temp)
#     stat_pop <- rbind(stat_pop, stat_pop_temp) 
#     
#   }
#   
#   
# } 
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
# stat_bloc <- stat_bloc %>%
#   mutate(P1.2 = case_when(P1=="WE"|P2=="WE" ~ "WE",
#                           (P1!= "WE" & P2 !="WE") & (P1=="EE" | P2 =="EE") ~ "EE",
#                           (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE") & (P1=="WA"|P2=="WA") ~ "WA",
#                           (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA") & (P2=="EA"|P2=="EA")~ "EA")) %>%
#   mutate(P2.2 = ifelse(P1==P1.2, as.character(P2), as.character(P1))) %>%
#   dplyr::select(-one_of("P1","P2")) %>%
#   rename(P1=P1.2, P2=P2.2) %>%
#   mutate(methode="haplotype") %>%
#   dplyr::select(bloc, P1, P2, FST, pos, nb_alleles, methode)
# 
# 
# write.table(stat_bloc, titre_stat_bloc, col.names=T, row.names=F, dec=".", sep="\t", quote=F)
# write.table(stat_pop, titre_stat_pop, col.names=T, row.names=F, dec=".", sep="\t", quote=F)
# 
# 
# 





sessionInfo()


