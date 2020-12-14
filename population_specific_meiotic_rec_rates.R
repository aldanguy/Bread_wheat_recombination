

# Extract population specific meiotic maps from PHASE outputs. 
# The hypothesis is that rho is proportionnal to c (here represented by Bayesian meiotic recombination rate of CsRe) and that the ratio between mean rho per genomic region (1AR1...7DR3) and and mean c per genomic region gives the proportionnality coefficient

Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))




cat("\n population_specific_meiotic_rec_rates.R \n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

titre_phase <- variables[1]
titre_csre <- variables[2]
titre_population_specific_meiotic_rec_rates <- variables[3]
titre_graphe <- variables[4]

# titre_phase <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/sorties_resumees_PHASE.txt"
# titre_csre <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_genetic_map.txt"
# titre_graphe <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/graphes/proportionnality_constant.png"


cat("\n\n Input 1 : Summary outputs from PHASE \n\n")
phase <- fread(titre_phase) 
head(phase)

phase <- phase %>%
  filter(w_center==T) %>% # removing windows border susceptible to have low quality inferences
  na.omit() %>% # removing intervals where PHASE failed to properly estimate historical recombination rates
  group_by(population, chr, region) %>%
  mutate(l=lpop*1e6) %>% # physical size of interval (in pb)
  mutate(d_historical=lambda_rho_med*l) %>% # historical genetic distance of interval. lambda_rho_med stands for "median of the posterior distribution of lambda*rho in this interval"
  summarise(d_historical=sum(d_historical), l= sum(l)) %>% # total genetic distance and total physical distance
  mutate(mean_lambda_rho_genomic_region=d_historical/l) %>%
  dplyr::select(population, chr, region, mean_lambda_rho_genomic_region)




cat("\n\n Input 2 : CsRe Bayesian meiotic recombination rates \n\n")
csre <- fread(titre_csre)
head(csre)

csre <- csre %>% group_by(chr, region) %>%
  mutate(l=lpop*1e6) %>%
  summarise(d_csre=sum(dbay), l=sum(l), mean_csre=d_csre/l) %>% # mean_csre is the average Bayesian meiotic recombination rate of CsRe per genomic region (in Morgan/pb)
  dplyr::select(chr, region, mean_csre)


coef_proportionnality <- phase %>% inner_join(csre, by=c("chr","region")) %>%
  mutate(coef_proportionnality=mean_lambda_rho_genomic_region/mean_csre) %>%
  dplyr::select(population, chr, region, coef_proportionnality)


population_specific <- fread(titre_phase)  %>%
  filter(w_center==T) %>% 
  na.omit() %>%
  dplyr::select(population, chr, region, posSNPlpop, posSNPrpop, SNPlpop, SNPrpop, lambda_rho_med) %>%
  inner_join(coef_proportionnality, by=c("population","chr","region")) %>%
  mutate(pop_specific_meiotic_rec_rate=lambda_rho_med/coef_proportionnality) %>%
  rename(posl=posSNPlpop, posr=posSNPrpop, SNPl=SNPlpop, SNPr=SNPrpop, historical_rec_rate=lambda_rho_med) %>%
  arrange(population, chr, posl) %>%
  mutate(meiotic_rec_rate=pop_specific_meiotic_rec_rate*1e8) %>% # conversion in cM/Mb
  dplyr::select(population, chr, region, posl, posr, SNPl, SNPr, historical_rec_rate, meiotic_rec_rate, coef_proportionnality)

# population_specific %>% group_by(population, chr) %>%
#   mutate(d_meiotic=pop_specific_meiotic_rec_rate*((posr-posl)/1e6)) %>%
#   summarise(d_meiotic=sum(d_meiotic)) # total genetic distance per chromosome (cM)

g <- coef_proportionnality %>%
  ungroup() %>%
  mutate(region=factor(region, levels = c("R1","R2a","C","R2b","R3"))) %>%
  mutate(population=factor(population, levels = c("WE","EE","WA","EA"))) %>%
  ggplot(aes(x=region, y=coef_proportionnality, col=population)) +
  geom_point() +
  facet_wrap(chr~., ncol=3, scales = "free") +
  xlab("") +
  ylab("Porportionnality constant between historical and CsRe meiotic rec rates") +
  scale_colour_manual(values=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA")) +
  theme_light() +
  stat_summary(fun.y=mean, colour="black", geom="line",group=1)
cat("\n\n Graphic for proportionnality constant \n\n")

ggsave(titre_graphe,g, width = 20, height = 20, units = "cm")

#ABDDA4 = green
#D7191C = blue
#FDAE61 = red
#2B83BA = orange


cat("\n\n Population specific meiotic recombination rates \n\n")
head(population_specific)
write.table(population_specific, titre_population_specific_meiotic_rec_rates, col.names = T, row.names = F, dec=".", sep="\t", quote=F)





sessionInfo()