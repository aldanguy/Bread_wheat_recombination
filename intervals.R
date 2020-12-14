Sys.time()
cat("\n\nintervals.R\n\n")
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
titre_pos <- variables[2]
titre_sortie <- variables[3]


# titre_mqs <- "/work/adanguy/these/pipeline/030420/temp/mqs.txt"                
# titre_pos <- "/work/adanguy/these/pipeline/030420/sources/SNP_positions.txt"   
# titre_sortie <- "/work/adanguy/these/pipeline/030420/temp/intervalles_de_novo.txt"


# titre_mqs <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/mqs.txt"
# titre_pos <-"/home/adanguydesd/Documents/These_Alice/pipeline/sources/SNP_positions.txt"
# titre_csre <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/mef/csre_mef.txt"


cat("\n\n Input 1 : physical positions of SNP \n\n")
pos <- fread(titre_pos, header=T, sep="\t", dec=".", data.table = F)
head(pos)

cat("\n\n Input 2 : polymorphic SNP in at least 1 population \n\n")
mqs <- fread(titre_mqs, header=F, sep="\t", dec=".", data.table = F) 
head(mqs)
mqs <- mqs %>%
  rename(SNP=V1) %>%
  unique() %>%
  left_join(pos, by="SNP") %>%
  arrange(chr, posSNP) %>%
  rename(SNPlint=SNP, posSNPlint=posSNP) %>%
  mutate(SNPrint=lead(SNPlint)) %>%
  mutate(posSNPrint=lead(posSNPlint)) %>%
  mutate(chrr=lead(chr)) %>%
  filter(chr == chrr) %>%
  dplyr::select(-d_cumul) %>%
  na.omit() %>%
  mutate(intervalle=paste0(chr,"_",posSNPlint,"_", posSNPrint))%>%
  mutate(lint=posSNPrint-posSNPlint) %>%
  mutate(lint=lint/1e6) %>%
  dplyr::select(chr, region, intervalle, SNPlint, SNPrint, posSNPlint, posSNPrint, lint)


cat("\n\n Output 1 : de-novo intervals \n\n")
head(mqs)
write.table(mqs, titre_sortie, col.names = T, row.names = F, dec=".", sep="\t", quote = F)


sessionInfo()