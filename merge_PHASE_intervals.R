Sys.time()
cat("\n\nmerge_PHASE_intervals.R\n\n")
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

titre_entree <- variables[1]
titre_sortie <- variables[2]
# 
# titre_entree <-  "/work/adanguy/these/pipeline/030420/fusion/echant_de_novo_intervals.txt"         
# titre_sortie <- "/work/adanguy/these/pipeline/030420/fusion/echant_de_novo_intervals_balanced.txt"

cat("\n\n Input 1 : Recombination estimates for the 5 populations within de-novo intervals \n\n")
data <- fread(titre_entree, header=T, sep="\t", dec=".", data.table = F)
head(data)

int <- table(data$intervalle, data$population) %>%
  as.data.frame() %>%
  spread(Var2,Freq) %>%
  filter(WA==1)%>%
  filter(EA==1)%>%
  filter(EE==1) %>%
  filter(WE==1) %>%
  dplyr::select(Var1) %>% 
  unlist() %>%
  as.vector()






data <- data %>% filter(intervalle %in% !! int) %>%
  arrange(population, chr, posSNPlint) 

data %>% group_by(chr, region) %>%
  filter(region != "C") %>%
  summarise(n=n()) %>%
  ungroup() %>%
  summarise(nmin=min(n), nmax=max(n))

cat("\n\n Output 1 : Recombination estimates for the 5 populations within de-novo intervals, but only-novo intervals with one estimate per population were kept \n\n")
head(data)
write.table(data, titre_sortie, col.names = T, row.names = F, dec=".", sep="\t", quote=F)

sessionInfo()
