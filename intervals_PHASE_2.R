
Sys.time()
cat("\n\nintervals_PHASE_2.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))





variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_mqs <- variables[1]
titre_phase <- variables[2]
titre_sortie <- variables[3]


# titre_mqs <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/int_polym.txt"
# titre_phase <-"/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/WE_1A_10_mef.txt"
# nb_echant <- 50
# nb_post <- 1000
# titre_sortie <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/echant/WE_3A_1_echant.txt"

# titre_mqs <- "/work/adanguy/these/pipeline/030420/temp/intervalles_de_novo.txt"                      
# titre_phase <- "/work/adanguy/these/pipeline/030420/temp/mef/EA_1A_9_sorties_brutes_PHASE.txt"       
# titre_sortie <- "/work/adanguy/these/pipeline/030420/temp/echant/WE_7D_108_echant_de_novo_intervals.txt"


# titre_mqs <- "/work/adanguy/these/pipeline/160320/temp/int_polym.txt"              
# titre_phase <- "/work/adanguy/these/pipeline/160320/temp/mef/WE_5A_153_mef.txt"      
# titre_sortie <- "/work/adanguy/these/pipeline/160320/temp/echant/WE_5A_153_echant.txt"


cat("\n\n Input 1 : Historical recombination estimates for the window \n\n")
phase <- fread(titre_phase, header=T, dec=".", sep="\t", data.table = F) 
head(phase)
phase <- phase %>%
  filter(w_center==TRUE) %>%
  filter(!is.na(PAC_slope)) %>%
  filter(!is.na(lambda_rho) & !is.na(rho) & !is.na(lambda)) %>%
  group_by(methode, population, chr, posSNPlpop, posSNPrpop, lpop, different_scaffold, SNPlpop, SNPrpop, ID) %>%
  summarize(lambda_rho=median(lambda_rho), lambda=median(lambda))

if (nrow(phase) > 0){
  
  cat("\n\n Input 2 : de-novo intervals \n\n")
  mqs <- fread(titre_mqs, header=T, dec=".", sep="\t", data.table = F) 
  print(head(mqs))
  mqs <- mqs %>%
    full_join(phase, by="chr") %>%
    filter(posSNPlint >= posSNPlpop & posSNPrint <= posSNPrpop) %>%
    arrange(posSNPlint) %>%
    rename(recombinaison=lambda_rho) %>%
    dplyr::select(population,
                  chr,
                  region,
                  ID,
                  intervalle,
                  SNPlint,
                  SNPlpop,
                  SNPrint,
                  SNPrpop,
                  posSNPlint,
                  posSNPlpop,
                  posSNPrint,
                  posSNPrpop,
                  lint,
                  lpop,
                  different_scaffold,
                  recombinaison,
                  lambda,
                  methode)
  
  cat("\n\n Output 1 : sampling in de-novo intervals \n\n")
  print(head(mqs))
  write.table(mqs, titre_sortie, col.names = T, row.names = F, dec=".", sep="\t", quote = F)
  
} else {
  
  cat("PAC_slope = NA")
  
}

sessionInfo()