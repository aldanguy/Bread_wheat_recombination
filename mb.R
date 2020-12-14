



Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))




cat("\n mb.R \n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

step <-  as.numeric(variables[1])
size <- as.numeric(variables[2])
titre_sortie <- variables[3]


# step <- 4e6 # parameters of resolution
# size <- 4e6


# titre_chr_tab <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Decoupage_chr_ble.tab"
# titre_phase <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/EA_1A_1_sorties_brutes_PHASE.txt"

# titre_phase <-  "/work/adanguy/these/pipeline/030420/fusion/EA_1A_sorties_brutes_PHASE.txt"
# titre_chr_tab <- "/work/adanguy/these/pipeline/amont/Decoupage_chr_ble.tab"                 
# titre_sortie <- "/work/adanguy/these/pipeline/030420/temp/Mb/carte_Mb_landraces_EA_1A.txt" 




i=1
chrs=paste0(rep(1:7, each=3), c("A","B","D"))

carte1M <- data.frame()
for (chr in chrs){
  
  
  cartetemp <- data.frame(chr=chr,posl=seq(1, 1e9-size, step)) %>%
    mutate(posr=posl+size-1)
  
  carte1M <- rbind(carte1M, cartetemp)
  
  
}


write.table(carte1M, titre_sortie, quote=F, row.names = F, dec=".", sep="\t")

sessionInfo()