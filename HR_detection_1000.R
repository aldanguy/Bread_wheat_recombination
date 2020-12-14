
Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bedr))




cat("\n HR_detection_1000.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_HR <- variables[1]
titre_HR_1000 <- variables[2]
# 
# titre_resume <-"/home/adanguydesd/Documents/These_Alice/pipeline/160320/sorties_resumees_PHASE.txt"
# titre_HR_100 <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/100_HR.txt"
# titre_HR <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/HR.txt"
# 
#   



cat("\n\n Input 1 : HR \n\n")
HR <- fread(titre_HR)

nb <- 1000


############  



# res2 <- HR %>%
#   filter(kept==T) %>%
#   group_by(population, chr) %>%
#   mutate(nHR=n()) %>%
#   mutate(ordre=seq(1,n(),1)) %>%
#   group_by(chr) %>%
#   mutate(nHR2=min(nHR)) %>%
#   ungroup() %>%
#   arrange(desc(lambda_med)) %>%
#   mutate(kept=ifelse(ordre <= nHR2, T, F)) %>%
#   filter(kept==T) %>%
#   group_by(population, chr, region) %>%
#   arrange(population, chr, posl) %>%
#   ungroup() %>%
#   dplyr::select(population, chr, region, posl, posr, taille, IDintHR, kept, lambda_med, d_genet_histo)
# 
res2 <- HR %>%
  filter(kept==T) %>%
  group_by(population) %>%
  mutate(nHR=n()) %>%
  ungroup() %>%
  mutate(nHR=min(nHR)) %>%
  group_by(population) %>%
  arrange(desc(lambda_med)) %>%
  mutate(kept=ifelse(row_number() <= min(nHR, nb), T, F)) %>%
  ungroup() %>%
  filter(kept==T) %>%
  arrange(population, chr, posl) %>%
  dplyr::select(population, chr, region, posl, posr, taille, IDintHR, kept, lambda_med, d_genet_histo)


table(res2$chr, res2$population)
table(res2$population)


cat("\n\n Output 1: 1000 most recombining interval per pop \n\n")
head(res2)
write.table(res2, titre_HR_1000, col.names = T, row.names = F, dec=".", sep="\t", quote=F)



sessionInfo()

