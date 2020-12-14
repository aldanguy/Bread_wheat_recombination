
Sys.time()
cat("\n\nmerge_PHASE_time.R\n\n")
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


titre_f <- variables[1]
titre_d <- variables[2]


 # titre_f <- "/work/adanguy/these/pipeline/170120/fusion/fenetres.txt"
 # titre_d <- "/work/adanguy/these/pipeline/170120/temp/autres/durees.txt"
 # titre_sortie <- titre_f

cat("\n\n Intput 1 : information about PHASE run for each window \n\n")
d <- fread(titre_d, header=F, dec=".", sep="\t", data.table = F)
head(d)
d <- d %>%
  rename(ID=V1, duree=V2, nb_underflows=V3, nb_fails=V4, PAC_slope=V5)

cat("\n\n Intput 2 : information about windows for PHASE \n\n")
f <- fread(titre_f, header=T, dec=".", sep="\t", data.table = F) 
head(f)
f <- f %>%
  dplyr::select(population,chr,region, poslfen,posrfen,centerlfen,centerrfen,nb_mqs,ID) %>%
  full_join(d, by="ID") %>%
  arrange(population, chr,poslfen)

cat("\n\n Output 1 : gathering information about windows \n\n")
head(f)
write.table(f, titre_f,col.names = T, dec=".",sep="\t", row.names = F, quote = F)
  


sessionInfo()