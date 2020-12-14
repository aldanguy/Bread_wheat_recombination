



### Extraction of polymorphic SNP in the 4 populations of landraces


Sys.time()
cat("\n\ncommon_SNP_2.R\n\n")
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

titre_WE_map <-variables[1]
titre_EE_map <-variables[2]
titre_WA_map <-variables[3]
titre_EA_map <-variables[4]
titre_sortie <- variables[5]


# titre_WE_map <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/WE.map"
# titre_EA_map <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/EA.map"

cat("\n\n Input 1 : PLINK file for WE polymorphic SNP \n")
WE <- fread(titre_WE_map)
head(WE)

cat("\n\n Input 2 : PLINK file for EE polymorphic SNP \n")
EE <- fread(titre_EE_map)
head(EE)

cat("\n\n Input 3 : PLINK file for WA polymorphic SNP \n")
WA <- fread(titre_WA_map)
head(WA)

cat("\n\n Input 4 : PLINK file for EA polymorphic SNP \n")
EA <- fread(titre_EA_map)
head(EA)


fusion <- rbind(WE,EE, WA, EA) %>%
  rename(chr=V1, SNP=V2, gen=V3, phy=V4) %>%
  group_by(SNP) %>%
  mutate(n=n()) %>%
  filter(n==4) %>%
  unique() %>%
  arrange(chr,phy) %>%
  dplyr::select(SNP) %>%
  unlist() %>%
  as.vector()

cat("\n\n Output 1 : PLINK file of polymorphic SNP in all populations \n")
head(fusion)
write.table(fusion, titre_sortie, col.names = F, row.names = F, dec=".", sep="\t", quote=F)

sessionInfo()


