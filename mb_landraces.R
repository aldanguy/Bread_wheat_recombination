



Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))




cat("\n mb_landraces.R \n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

titre_phase <- variables[1]
titre_chr_tab <- variables[2]
titre_fenetres <- variables[3]
titre_sortie <- variables[4]


 # step <- 4e6  parameters of resolution
 # size <- 4e6


 # titre_chr_tab <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Decoupage_chr_ble.tab"
 # titre_phase <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/EA_1A_1_sorties_brutes_PHASE.txt"
 # 
# titre_phase <-  "/work/adanguy/these/pipeline/030420/fusion/WA_2B_sorties_brutes_PHASE.txt"
# titre_chr_tab <- "/work/adanguy/these/pipeline/amont/Decoupage_chr_ble.tab"                 
# #  # titre_sortie <- "/work/adanguy/these/pipeline/030420/temp/Mb/carte_Mb_landraces_EA_1A.txt" 
# # 
# # titre_phase <- "/work/adanguy/these/pipeline/030420/fusion/EA_1A_sorties_brutes_PHASE.txt"      
# titre_fenetres <- "/work/adanguy/these/pipeline/030420/temp/carte_chrregion.txt"                   
# # titre_sortie <-  "/work/adanguy/these/pipeline/030420/temp/Mb/carte_chrregion_landraces_EA_1A.txt"
# # 


cat("\n\n Input 1 : Raw outputs from PHASE \n\n")
phase <- fread(titre_phase)
head(phase)
phase <- phase %>%
  filter(w_center==T) %>%
  na.omit() %>%
  arrange(chr, posSNPlpop, num_row_PHASE_output) %>%
  mutate(lpop=lpop*1e6) %>%
  dplyr::select(population, chr, ID, posSNPlpop, posSNPrpop, lpop, lambda_rho, num_row_PHASE_output)




cat("\n\n Input 2 : windows \n\n")
carte1M <- fread(titre_fenetres)
head(carte1M)

regions  <- read.table(titre_chr_tab, header=T, dec=".", sep="\t", skip=1) %>%
  dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
  mutate(chr = str_remove(Chromosome,"chr")) %>%
  dplyr::select(-one_of("Chromosome"))

minmax <- phase %>%
  group_by(chr) %>% summarise(min=min(posSNPlpop), max=max(posSNPrpop))


carte1M <- carte1M %>% inner_join(minmax, by="chr") %>%
  filter(!(posl >= max | posr <= min)) %>%
  dplyr::select(chr, posl, posr) %>%
  arrange(chr, posl)



i=1
k=0
for (i in 1:nrow(carte1M)){
  
  
  phase2 <- phase %>% filter(!(posSNPrpop <= !!carte1M$posl[i] | posSNPlpop >= !!carte1M$posr[i] ))
  
  if(nrow(phase2)>0){
    k=k+1
    
     # phasetemp <- carte1M %>%
     #   slice(i) %>%
     #   full_join(phase2, by="chr") %>%
     #   filter(!(posSNPlpop >=posr | posSNPrpop <= posl)) %>%
     #   group_by(population, chr, posl, posr, ID, num_row_PHASE_output, posSNPlpop) %>%
     #   mutate(L1=min(posSNPrpop, posr) - max(posSNPlpop, posl) ) %>% 
     #   mutate(d=lambda_rho*L1) %>%
     #   group_by(population, chr,posl, posr, ID, num_row_PHASE_output) %>%
     #   summarise(d2=sum(d), l2=sum(L1), lambda_rho=d2/l2) %>%
     #   group_by(population, chr, posl, posr, ID, l2) %>%
     #   summarise(lambda_rho=mean(lambda_rho)) %>%
     #   mutate(d3=lambda_rho*l2) %>%
     #   group_by(population, chr, posl, posr) %>%
     #   summarise(d4=sum(d3), l4=sum(l2), lambda_rho=d4/l4) %>%
     #   full_join(regions, by="chr") %>%
     #   mutate(region = ifelse(posl <=  R1.R2a*1e6, "R1",
     #                          ifelse(posl >= R1.R2a*1e6 +1 & posl <= R2a.C*1e6, "R2a",
     #                                 ifelse(posl >= R2a.C*1e6 +1 & posl <= C.R2b*1e6, "C", 
     #                                        ifelse(posl >= C.R2b*1e6 +1 & posl <= R2b.R3*1e6, "R2b",
     #                                               ifelse(posl >= R2b.R3*1e6 +1, "R3", NA)))))) %>%
     #   arrange(chr, posl) %>%
     #   na.omit() %>%
     #   dplyr::select(population, chr, region, posl, posr, lambda_rho) %>%
     #   as.data.frame()
     # 
    
    phasetemp <- carte1M %>%
      slice(i) %>%
      full_join(phase2, by="chr") %>%
      filter(!(posSNPlpop >=posr | posSNPrpop <= posl)) %>%
      group_by(population, chr, posl, posr, ID, num_row_PHASE_output, posSNPlpop) %>%
      mutate(L1=min(posSNPrpop, posr) - max(posSNPlpop, posl) ) %>% 
      mutate(d=lambda_rho*L1) %>%
      group_by(population, chr,posl, posr, ID, num_row_PHASE_output) %>%
      summarise(d2=sum(d), l2=sum(L1), lambda_rho=d2/l2) %>%
      group_by(population, chr, posl, posr, num_row_PHASE_output) %>%
      summarise(d3=sum(d2), l3=sum(l2), lambda_rho=d3/l3) %>%
      group_by(population, chr, posl, posr) %>%
      summarise(sd_lambda_rho=sd(log10(lambda_rho)), lambda_rho=mean(lambda_rho)) %>%  # ou sd_lambda_rho=sd(lambda_rho)
      full_join(regions, by="chr") %>%
      mutate(region = ifelse(posl <=  R1.R2a*1e6, "R1",
                             ifelse(posl >= R1.R2a*1e6 +1 & posl <= R2a.C*1e6, "R2a",
                                    ifelse(posl >= R2a.C*1e6 +1 & posl <= C.R2b*1e6, "C", 
                                           ifelse(posl >= C.R2b*1e6 +1 & posl <= R2b.R3*1e6, "R2b",
                                                  ifelse(posl >= R2b.R3*1e6 +1, "R3", NA)))))) %>%
      arrange(chr, posl) %>%
      na.omit() %>%
      dplyr::select(population, chr, region, posl, posr, lambda_rho, sd_lambda_rho) %>%
      as.data.frame()
    
    
    if (k==1){
      
      
      
      cat("\n\n Landraces historical map with specific resolution\n\n")
      print(head(phasetemp) )
      write.table(phasetemp, titre_sortie, col.names = T, row.names = F, dec=".", sep="\t", quote=F, append = F)
      
      
      
      
    } else {
      
      
      write.table(phasetemp, titre_sortie, col.names = F, row.names = F, dec=".", sep="\t", quote=F, append = T)
      
      
      
    }
    
    
  }  
}


sessionInfo()