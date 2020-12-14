Sys.time()
cat("\n\nlandraces.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(janitor))






variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")



titre_l <-variables[1]
titre_s <-variables[2]
titre_sortie <- variables[3]
titre_K4 <-variables[4]
titre_K8 <-variables[5]

  
      # titre_l <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/accessions_4403_alice.csv"
      # titre_s <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/structure_tot_def.txt"
      # titre_sortie <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/landraces.txt"
      # titre_K4 <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/autre_amont/structure/landrace8741hap_k4_rep5_f_f-ord-ord-ord"
      # titre_K8 <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/autre_amont/structure/landrace8741hap_k8_rep2_f_f-ord-ord-ord"



country <- read.table(titre_l, header=T, dec=".", sep=";", colClasses = "character") %>%
  dplyr::select(country_code, country) %>%
  unique() %>%
  arrange(country_code) %>%
  na.omit()


l <-  read.table(titre_l, header=T, dec=".", sep=";", colClasses = "character") %>%
  filter(landrace==1)

s <- read.table(titre_s, header=T, dec=".", sep="\t") %>%
  mutate(origine_genet=as.character(origine_genet)) %>%
  mutate(ERGE = as.character(ERGE)) %>%
  dplyr::select(Unit,Nom_TaBW420K, ERGE, LINE, origine_genet, memb4, gp4, memb8, gp8, region, color) %>%
  rename(country_code=origine_genet, STRmemb4_max=memb4, STRgp4=gp4, STRmemb8_max=memb8, STRgp8=gp8, STRarea8=region, STRcolor8=color)%>%
  left_join(country, by="country_code") %>%
  full_join(l, by="ERGE", suffix=c("",".x")) %>%
  dplyr::select(-LINE) %>%
  rename(LINE=LINE.x)%>%
  dplyr::select(Unit,Nom_TaBW420K, ERGE, LINE, country, STRmemb4_max, STRgp4, STRmemb8_max, STRgp8, STRarea8, STRcolor8) %>%
  arrange(Unit) 


k4 <- read.table(titre_K4, header=F, dec=".", sep="\t") %>%
  remove_constant() 

colnames(k4) <- paste0("STRmemb4_", 1:4)
  

k8 <- read.table(titre_K8, header=F, dec=".", sep="\t") %>%
  remove_constant() 

colnames(k8) <- paste0("STRmemb8_", 1:8)

s <- cbind(s,k4,k8) 


head(s)

# From an other passeport data file
s[which(s$Nom_TaBW420K=="86_ERGE24052_F11_121213"),"country"] <- "Italy"
  
write.table(s,titre_sortie, col.names = T, dec=".",sep="\t", row.names = F, quote = F)


sessionInfo()

