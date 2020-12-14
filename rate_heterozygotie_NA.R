


Sys.time()
cat("\n\nrate_heterozygotie_NA.R\n\n")
rm(list = ls())
set.seed(1)
graphics.off()


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(trio))






variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


pop <- variables[1]
titre_ped <- variables[2]
titre_map <- variables[3]
titre_freq <- variables[4]
titre_pos <- variables[5]
titre_SNP_PLINK <- variables[6]
titre_landraces_to_keep <- variables[7]
titre_stat_locus <- variables[8]
titre_stat_ind <- variables[9]
seuil_het <- as.numeric(variables[10])
seuil_NA <- as.numeric(variables[11])
seuil_maf <-  as.numeric(variables[12])

  # pop <- "landraces"         
  #  titre_map <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/Map.map"
  #  titre_ped <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/Ped.ped"
  #  titre_landraces_to_keep <- "/work/adanguy/these/pipeline/140220/temp/autres/landraces_to_keep.txt"
  #  titre_pos <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/SNP_positions.txt"
  #  seuil_het <- 0.05
  #  seuil_NA <- 0.1
  #  seuil_maf <- 1
  

# pop <-  "landraces"                                                            
# titre_ped <- "/work/adanguy/these/pipeline/030420/sources/landraces.ped"            
# titre_map <- "/work/adanguy/these/pipeline/030420/sources/landraces.map"            
# titre_freq <- "/work/adanguy/these/pipeline/030420/sources/landraces.frq"            
# titre_pos <- "/work/adanguy/these/pipeline/030420/sources/SNP_positions.txt"        
# titre_SNP_PLINK <- "/work/adanguy/these/pipeline/030420/temp/autres/SNP_PLINK_2.txt"      
# titre_landraces_to_keep <-"/work/adanguy/these/pipeline/030420/temp/autres/landraces_to_keep.txt"
# titre_stat_locus <-"/work/adanguy/these/pipeline/030420/temp/autres/NA.txt"               
# titre_stat_ind <- "/work/adanguy/these/pipeline/030420/temp/autres/NA.txt"               
# seuil_het <-  "0.05"                                                                 
# seuil_NA <- "1"                                                                    
# seuil_maf <- "0" 




cat("\n\n Input 1 : Physical position of SNP \n\n")
pos <- fread(titre_pos , header=T, dec=".", sep="\t", data.table = F)
head(pos)

cat("\n\n Input 2 : SNP of landraces data \n\n")
map <- fread(titre_map, header=F, dec=".", sep="\t", data.table = F) 
head(map)
map <- map %>%
  rename(chr_code_nombre=V1, SNPpop=V2, pos_gen=V3, posSNP=V4) %>%
  dplyr::select(SNPpop)

cat("\n\n Input 3 : frequency SNP of landraces data \n\n")
d <- fread(titre_freq, header=T, sep=" ", dec=".", data.table = F) 
head(d)
d <- d %>%
  dplyr::select(SNP, MAF) %>%
  rename(SNPpop=SNP)


# Select 50 individuals and 100 locus
# bilan <- read.pedfile(titre_ped, first.row = T) %>% dplyr::select(-one_of("famid","fatid","motid","sex","affected", "pid")) %>% slice(1:50)
# bilan <- bilan[,c(1:200)]
# map <- fread(titre_map, header=F, dec=".", sep="\t", data.table = F) %>% rename(chr_code_nombre=V1, SNPpop=V2, pos_gen=V3, posSNP=V4) %>% slice(1:200)
# 

cat("\n\n Input 4 : genotyped SNP of landraces data \n\n")
bilan <- read.pedfile(titre_ped, first.row = T) 
bilan[1:10,1:10]

noms_individus <- bilan %>% dplyr::select(pid) %>%
  unlist() %>%
  as.vector()

bilan2 <- bilan %>%
  dplyr::select(-one_of("famid","fatid","motid","sex","affected", "pid")) %>%
  t() %>%
  as.data.frame() %>%
  mutate_at(vars(everything()), funs(str_replace_all(., pattern = c("0"), replacement = "N"))) %>%
  mutate_at(vars(everything()), funs(paste0(., lead(.)))) %>%
  slice(seq(1,n(),2)) %>%
  mutate_at(vars(everything()), funs(str_replace_all(., pattern = c("AA|TT|GG|CC"), replacement = "0"))) %>%
  mutate_at(vars(everything()), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'T','G','C')), 2))), collapse = "|"), replacement= "1"))) %>%
  mutate_at(vars(everything()), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'T','G','C','N')), 2))), collapse = "|"), replacement= "NA"))) %>%
  mutate_at(vars(everything()), funs(as.numeric(as.character(.)))) 


locus2 <- bilan2 %>%
  as.data.frame() %>%
  mutate(n = rowSums(.[grep("[0-9]", names(.))], na.rm = TRUE)) %>%
  mutate(L = rowSums(!is.na(.[grep("[0-9]", names(.))]))) %>%
  mutate(taux_het = n/L) %>%
  dplyr::select(taux_het) %>%
  cbind(., map) %>%
  left_join(pos, by=c("SNPpop"="SNP")) %>%
  left_join(d, by="SNPpop") %>%
  mutate(population=!!pop) %>%
  dplyr::select(population, chr,chr_code_nombre, region, SNPpop, posSNP,other_SNP_ID,category, taux_het, MAF) 
# dplyr::select(population, chr,chr_code_nombre, region, SNPpop, posSNP,other_SNP_ID,category, taux_het) 




individus2 <- bilan2 %>%
  t() %>%
  as.data.frame() %>%
  mutate(n = rowSums(.[grep("V", names(.))], na.rm = TRUE)) %>%
  mutate(L = rowSums(!is.na(.[grep("V", names(.))]))) %>%
  mutate(taux_het = n/L) %>%
  mutate(population=!!pop) %>%
  mutate(LINE=!!noms_individus) %>%
  dplyr::select(LINE, taux_het)





rm(bilan2)
rm(map)
rm(d)
rm(pos)





bilan3 <- bilan %>%
  dplyr::select(-one_of("famid","fatid","motid","sex","affected", "pid")) %>%
  t() %>%
  as.data.frame() %>%
  mutate_at(vars(everything()), funs(str_replace_all(., pattern = c("A|T|G|C"), replacement = "1"))) %>%
  mutate_at(vars(everything()), funs(as.numeric(as.character(.)))) 



locus3 <- bilan3 %>%
  as.data.frame() %>%
  mutate(n = rowSums(.[grep("V", names(.))], na.rm = TRUE)) %>%
  mutate(ncol=length(.[grep("V", names(.))])) %>%
  dplyr::select(n, ncol) %>%
  mutate(locus=as.factor(rep(seq(1,nrow(locus2),1), each=2))) %>%
  group_by(locus) %>%
  summarise(taux_NA = 1-(sum(n)/sum(ncol)), test=sd(n)) %>%
  dplyr::select(taux_NA) %>%
  cbind(locus2,.)

 
individus3 <- bilan3 %>%
  t() %>%
  as.data.frame() %>%
  mutate(n = rowSums(.[grep("V", names(.))], na.rm = TRUE)) %>%
  mutate(ncol=length(.[grep("V", names(.))])) %>%
  dplyr::select(n, ncol) %>%
  mutate(taux_NA = 1 - (n/ncol)) %>%
  dplyr::select(taux_NA) %>%
  cbind(individus2,.)



SNP_to_keep <- locus3 %>%
  filter(taux_het <= seuil_het & taux_NA <= seuil_NA & MAF >=seuil_maf) %>%
  dplyr::select(SNPpop) %>%
  unlist() %>%
  as.vector()




landraces_to_keep <- individus3 %>%
  filter(taux_het <= seuil_het & taux_NA <= seuil_NA) %>%
  dplyr::select(LINE) %>%
  arrange(LINE) %>%
  unlist() %>%
  as.vector()


cat("\n\n Output 1 : % het, % missing data of landraces \n\n")
write.table(individus3, titre_stat_ind, sep="\t", col.names = T, row.names = F, dec=".", quote=F)
head(individus3)

cat("\n\n Output 2 : % het, % missing data, MAF of SNP \n\n")
write.table(locus3, titre_stat_locus, sep="\t", col.names = T, row.names = F, dec=".", quote=F)
head(locus3)

cat("\n\n Output 3 :list of SNP to keep after filtering \n\n")
write.table(SNP_to_keep, titre_SNP_PLINK, sep="\t", col.names = F, row.names = F, dec=".", quote=F)
head(SNP_to_keep)

cat("\n\n Output 4 :list of landraces to keep after filtering \n\n")
write.table(landraces_to_keep, titre_landraces_to_keep, sep="\t", col.names = F, row.names = F, dec=".", quote=F)
head(landraces_to_keep)


sessionInfo()

