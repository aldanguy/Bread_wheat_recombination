
#### Formatting of tables

Sys.time()
cat("\n\nnumeric_results.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(cocor))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(nlme))
suppressPackageStartupMessages(library(weights))



variables <- commandArgs(trailingOnly=TRUE)

cat("\n Graphics and stats\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

titre_map4mb <- variables[1]
titre_hris <-variables[2]
titre_phase <- variables[3]
titre_cormm <- variables[4]
titre_cliques <- variables[5]
titre_cormm_common <- variables[6]
titre_fst <- variables[7]




titre_map4mb <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/historical_and_meiotic_maps_4Mb_published.txt"
titre_hris <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/HR_published.txt"
titre_phase <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/PHASE_summary_outputs_published.txt"
titre_cormm <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/supplementary_file_S8_correlations_mixed_models.txt"
titre_cliques <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/cliques_published.txt"
titre_cliques <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/cliques_published.txt"
titre_fst <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/FST_published.txt"





head(fread(titre_map4mb))
head(fread(titre_hris))
head(fread(titre_phase))
head(fread(titre_cormm))
head(fread(titre_cliques))
head(fread(titre_fst))

#### Zou's test
# The two correlations are overlapping, i.e., they have one variable in common. 
# This test calculates the confidence interval of the difference between the two correlation coefficients (ex: r_{WECsRe} and r_{EECSRE}). 
# Because the significance depends on the intercorrelation between WE and EE (r.WEEE), this intercorrelation has to be provided as an additional parameter. 
# If the confidence interval includes zero, the null hypothesis that the two correlations are equal must be retained. 
# If the confidence interval does not include zero, the null hypothesis has to be rejected.



map4mb <- fread(titre_map4mb) %>% 
  pivot_wider(id_cols = c("chr", "posl", "posr"), values_from = "rec_rate", names_from = "population") %>%
  na.omit() %>%
  as.data.frame()


head(map4mb)


corWEEECSRE <- cocor(~ WE + CsRe | EE + CsRe, data=map4mb, return.htest = T, test="zou2007", conf.level =1-1e-3)
corWEWACSRE <- cocor(~ WE + CsRe | WA + CsRe, data=map4mb, return.htest = T, test="zou2007", conf.level =1-1e-3)
corWEEACSRE <- cocor(~ WE + CsRe | EA + CsRe, data=map4mb, return.htest = T, test="zou2007", conf.level =1-1e-16)
corEEWACSRE <- cocor(~ EE + CsRe | WA + CsRe, data=map4mb, return.htest = T, test="zou2007", conf.level =1-0.4)
corEEEACSRE <- cocor(~ EE + CsRe | EA + CsRe, data=map4mb, return.htest = T, test="zou2007", conf.level =1-1e-10)
corWAEACSRE <- cocor(~ WA + CsRe | EA + CsRe, data=map4mb, return.htest = T, test="zou2007", conf.level =1-1e-9)

# confidence intervals
corWEEECSRE$zou2007$conf.int[1:2] # 1e-4 < p.value < 1e-3
corWEWACSRE$zou2007$conf.int[1:2] # 1e-4 < p.value > 1e-3
corWEEACSRE$zou2007$conf.int[1:2] # p.value < 1e-16
corEEWACSRE$zou2007$conf.int[1:2] # 0.4 < p.value < 0.5
corEEEACSRE$zou2007$conf.int[1:2] # 1e-11 < p-value < 1e-10
corWAEACSRE$zou2007$conf.int[1:2] # 1e-10 < p-value < 1e-9




### density of HRIs and CsRe recombination rate


csre4mb <- fread(titre_map4mb) %>%
  filter(region !="C") %>%
  filter(population=="CsRe")


hris <- fread(titre_hris) %>%
  dplyr::select(-lambda_med) %>%
  inner_join(csre4mb, by="chr", suffix=c(".landraces",".csre")) %>%
  filter(!(posr.csre <= posl.landraces | posr.landraces <= posl.csre)) %>%
  group_by(population.landraces, chr, region.csre, posl.csre, posr.csre, rec_rate) %>%
  summarise(nhris=length(unique(HRI_ID))) 

head(hris)

summary(lm(rec_rate ~ nhris + region.csre, data=hris))
# effect of number of HRIs per intervals have a >0 and significant effect on CsRe rec rate

###density of HRIs per region

phase <- fread(titre_phase) %>%
  filter(method="population_specific_SNPs") %>%
  filter(region !="C") %>%
  group_by(region) %>%
  summarise(nintervals=n())

phris <- fread(titre_hris) %>%
  group_by(region) %>%
  summarise(nhris=n()) %>%
  inner_join(phase, by=c("region")) %>%
  mutate(p = round(nhris/nintervals,2))
phris



chisq_test_0 <- chisq.test(x=phris %>% dplyr::select(-region, -p) %>% as.matrix())
chisq_test_0
chisq_test_0$expected

# density of HRIs in populations


phase <- fread(titre_phase) %>%
  filter(method="population_specific_SNPs") %>%
  filter(region !="C") %>%
  group_by(population) %>%
  summarise(nintervals=n())
head(phase)




phris <- fread(titre_hris) %>%
  group_by(population) %>%
  summarise(nhris=n()) %>%
  inner_join(phase, by=c("population")) %>%
  mutate(p = round(nhris/nintervals,2))


phris




# HRIs represent from 1% of intervals in EA to 2% in WE



chisq_test <- chisq.test(x=phris %>% dplyr::select(-population, -p) %>% as.matrix())
chisq_test
chisq_test$expected
# EA and EE have a deficit in HRIS while WE and WA has more than expected




### Differences in average correlation per pair of populations

cormm <- fread(titre_cormm) %>%
  filter(method=="population_specific_SNPs") %>%
  filter(region !="C") %>%
  mutate(pair_of_pop=paste0(P1,P2))

pairwise_t_test <- pairwise.t.test(x=cormm$cor_lambda, g=cormm$pair_of_pop, p.adjust.method = "bonf")
pairwise_t_test
matrice <- as.matrix(pairwise_t_test$p.value)
matrice=rbind.data.frame(NA,matrice)
matrice=cbind.data.frame(matrice,NA)
colnames(matrice) <- c(colnames(matrice)[1], row.names(matrice)[-1])
row.names(matrice) <- c(colnames(matrice)[1], row.names(matrice)[-1])
diag(matrice) <- 1
matrice <- as.matrix(matrice)
matrice[which(matrice < 0.05)] <- 0
matrice[which(matrice >= 0.05)] <- 1
plot(graph_from_adjacency_matrix(matrice, mode = "lower", weighted = T, diag = F),vertex.size=3, main="not significant differences")
# to put in relationship with
cormm %>% group_by(pair_of_pop) %>%
  summarise(mean=mean(cor_lambda)) %>%
  arrange(desc(mean))
min(matrix(pairwise_t_test$p.value, nrow=1), na.rm=T)
sort(matrix(pairwise_t_test$p.value, nrow=1))


#overall test
summary(lm(cor_lambda ~ pair_of_pop, data = cormm))


#### densiy of shared hot windows along genome does not depend on genomic region
dhris2 <- fread(titre_hris) %>%
  group_by(region) %>%
  summarise(nhris=n())

dcliques_shared <- fread(titre_cliques) %>%
  filter(!populations %in% c("WE","EE","WA","EA")) %>%
  group_by(region) %>%
  summarise(ncliques=n()) %>%
  inner_join(dhris2, by="region")
dcliques_shared

chisq_test_2 <- chisq.test(x=dcliques_shared %>% dplyr::select(-region) %>% as.matrix())
chisq_test_2
chisq_test_2$expected



# intensity in shared windows increased
intensity <- fread(titre_cliques) %>%
  mutate(shared=case_when(populations %in% c("WE","EE","WA","EA") ~ "1",
                          populations %in% c("WEEE","WEWA","WEEA","EEWA","EEEA","WAEA")~ "2",
                          populations %in% c("WEEEWA","WEEEEA","WEWAEA","EEWAEA" )~ "3",
                          populations %in% c("WEEEWAEA")~ "4")) %>%
  inner_join(fread(titre_hris), by="chr", suffix=c(".cliques",".hris")) %>%
  filter(!(posr.cliques < posl.hris | posr.hris < posr.cliques)) %>%
  group_by(shared) 
  
intensity %>%
  summarise(median=round(median(lambda_med),1))
pairwise_wilcox_test <- pairwise.wilcox.test(intensity$lambda_med, intensity$shared, p.adjust.method = "bonf")  
pairwise_wilcox_test
min(matrix(pairwise_wilcox_test$p.value, nrow=1), na.rm=T)
max(matrix(pairwise_wilcox_test$p.value, nrow=1), na.rm=T)




### decreasing relationship with correlation in pop.specific dataset and common SNP dataset


cormm_fst <- fread(titre_cormm) %>% 
  inner_join(fread(titre_fst) %>% filter(method=="haplotypic_blocks"), by= c("chr","region","P1","P2")) %>%
  filter(region !="C") %>%
  mutate(pair_of_pop=paste0(P1,P2)) %>%
  mutate(chrregion=paste0(chr,region))

cormm_common_fst <- fread(titre_cormm_common) %>% 
  inner_join(fread(titre_fst) %>% filter(method=="haplotypic_blocks"), by= c("chr","region","P1","P2")) %>%
  filter(region !="C") %>%
  mutate(pair_of_pop=paste0(P1,P2))%>%
  mutate(chrregion=paste0(chr,region))

mod = lmList( cor_lambda ~ FST | chrregion, data=cormm_fst)
mod_common = lmList( cor_lambda ~ FST | chrregion, data=cormm_common_fst)
sd <- as.numeric()
sd_common <- as.numeric()
effect <- as.numeric()
effect_common <- as.numeric()
for (chrregion in unique(cormm_fst$chrregion)){

  
effect <- c(effect, summary(mod[[chrregion]])$coefficients[2,1])
effect_common <- c(effect_common, summary(mod_common[[chrregion]])$coefficients[2,1])

sd <- c(sd, summary(mod[[chrregion]])$coefficients[2,2])
sd_common <- c(sd_common, summary(mod_common[[chrregion]])$coefficients[2,2])

  
}
wtd_t_test <- wtd.t.test(x=effect, y=effect_common, weight = sd, weighty = sd_common, bootse=TRUE)
wtd_t_test


sessionInfo()
