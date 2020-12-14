
Sys.time()
cat("\n\ncsre_genotyping.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)

suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))






variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

titre_matrice <-variables[1]
titre_autres_SNP <-variables[2]
titre_pos <-variables[3]
titre_output <-variables[4]


  # titre_matrice <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/MatriceGenoHead.CorrectedName.UniqSorted.tab"
  # titre_autres_SNP <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/CartPolyHigh-Oriented_def_def_def.txt"
  # titre_pos <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/SNP_positions.txt"
  # titre_output <- "/home/adanguydesd/Documents/These_Alice/pipeline/sources/genotypage_csre.txt"
  # 
  # 

  # titre_matrice <- "/work/adanguy/these/pipeline/amont/MatriceGenoHead.CorrectedName.UniqSorted.tab"
  # titre_autres_SNP <-"/work/adanguy/these/pipeline/amont/CartPolyHigh-Oriented_def_def_def.txt"       
  # titre_pos <- "/work/adanguy/these/pipeline/160320/sources/SNP_positions.txt"                  
  # titre_output <- "/work/adanguy/these/pipeline/160320/sources/matrice_genotypage_csre.txt"


pos <- fread(titre_pos, header=T, dec=".", sep="\t", data.table = F) %>%
  dplyr::select(one_of("SNP", "chr", "region", "posSNP", "other_SNP_ID"))





SNP_published <- fread(titre_autres_SNP, header=T, dec=".", sep="\t", data.table = F) %>%
  filter(!is.na(BW_published)) %>%
  dplyr::select(one_of("BW_published")) %>%
  unlist()%>%
  as.vector()




  
matrice <- fread(titre_matrice, header=F, sep="\t", dec=".") %>%
  `colnames<-`(as.vector(unlist(read.table(titre_matrice, sep=" ", dec=".", header=F, nrow=1)))) %>%
  rename(other_SNP_ID = Snp_name) %>%
  mutate(other_SNP_ID = gsub("\\*","",other_SNP_ID)) %>%
  left_join(pos, by="other_SNP_ID") %>%
  filter(SNP %in% !!SNP_published) %>%
  dplyr::select(-one_of("other_SNP_ID")) %>%
  mutate(population="CsRe") %>%
  dplyr::select(population, chr, region, SNP, posSNP, everything()) %>%
  arrange(chr, posSNP)
  
matrice[1:10,1:6]

  
write.table(matrice, titre_output, sep="\t", col.names = T, row.names = F, dec=".", quote=F)

sessionInfo()