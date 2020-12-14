


# Study the relationship between the average CsRe rec rate and the number of HR of landraces population within 4 Mb windows



Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))




cat("\n csre_and_HR.R \n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

titre_hr <- variables[1]
titre_csre <- variables[2]
titre_co <- variables[3]
titre_graphe <- variables[4]

# titre_csre<- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/map_4mb_csre.txt"
# titre_hr <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/HR.txt"
# titre_graphe <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/graphes/association_intensity_density.png"
# titre_co <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/csre_crossovers_positions.txt"


cat("\n\n Average CsRe recombination rate per 4 Mb window \n")
carte <- fread(titre_csre) 
head(carte)

cat("\n\n HR \n")
hr <- fread(titre_hr)
head(hr)


cat("\n\n csre crossovers positions \n")
co <- fread(titre_co)
head(co)

carte <- carte %>%
  na.omit() %>%
  rename(poslw=posl, posrw=posr) %>% 
  filter(region!="C")
  
hr <- hr %>% filter(kept==T) %>%
  inner_join(carte, by=c("chr","region")) %>%
  rename(population=population.x) %>%
  group_by(population, chr, region, poslw, posrw, cbay) %>%
  summarise(n=length(which(!(poslw >= posr | posrw <= posl)))) %>%
  ungroup() %>%
  mutate(population=factor(population, levels = c("WE","EE","WA","EA")))
head(hr)

tab <- data.frame()

for (r in c("R1","R2a","R2b","R3")) {
  
  mod <- lm(cbay~n+population+n*population, data=hr %>% filter(region==!!r))
  
  temp <- data.frame(region=r,
                     coefficient=mod$coefficient[2],
                     p_value=summary(mod)$coefficient[2,4],
                     diff_pop=min(summary(mod)$coefficient[6:8,4]) < 0.05)
  tab <- rbind(tab, temp)
  
  
}

tab


graphe <- ggplot(hr, aes(x = as.factor(n), y=cbay, col=population)) + geom_boxplot() +
  facet_grid(population~region, scales = "free") +
  geom_point() +
  scale_y_continuous(trans="log10") +
  xlab("Number of HR within 4 Mb window") +
  scale_colour_manual(values=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA"))+
  ylab("CsRe Bayesian recombination rate per 4 Mb window (cM/Mb)") +
  theme_light() 

graphe

cat("\n\n Graphe relationship between HR and CsRe rec rates \n")
ggsave(titre_graphe, graphe)



co <- fread(titre_co) %>% group_by(chr, IDco) %>% 
  summarise(posl=min(posSNPlpop), posr=max(posSNPrpop)) %>%
  ungroup() %>%
  mutate(nbCO=n())

# proportion of CsRe crossover falling within HR
fread(titre_hr) %>% filter(kept==T) %>%
  dplyr::select(population, chr, posl, posr, IDintHR) %>%
  full_join(co, by="chr") %>%
  filter(!(posl.x >= posr.y | posr.x <= posl.y)) %>%
  dplyr::select(population, IDco, nbCO) %>%
  unique() %>%
  group_by(population,nbCO) %>%
  summarise(pop_co_HR=unique(n()/nbCO)) %>%
  dplyr::select(population, pop_co_HR)

# proportion of HR falling within a CO
fread(titre_hr) %>% filter(kept==T) %>%
  group_by(population) %>%
  mutate(nbHR=n()) %>%
  ungroup() %>%
  dplyr::select(population, chr, posl, posr, IDintHR, nbHR) %>%
  full_join(co, by="chr") %>%
  filter(!(posl.x >= posr.y | posr.x <= posl.y)) %>%
  dplyr::select(population, IDintHR, nbHR) %>%
  unique() %>%
  group_by(population,nbHR) %>%
  summarise(pop_HR_co=unique(n()/nbHR)) %>%
  dplyr::select(population, pop_HR_co, nbHR)
  
  



sessionInfo()
