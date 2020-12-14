


Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MESS))
suppressPackageStartupMessages(library(ggpubr))




cat("\n gini_populations.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_resume <-variables[1]
titre_echant <- variables[2]
titre_csre <- variables[3]
titre_sortie_temp <- variables[4]
titre_graphe <- variables[5]
titre_sortie<- variables[6]


 # titre_graphe <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/graphes/gini.png"
 # titre_resume <-"/home/adanguydesd/Documents/These_Alice/pipeline/160320/sorties_resumees_PHASE.txt"
 # titre_echant <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/echant_de_novo_intervals_balanced.txt"
 # titre_sortie_temp <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/temporaire.txt"
 # titre_csre <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_genetic_map.txt"

# titre_resume <- "/work/adanguy/these/pipeline/020820/PHASE/PHASE_summary_outputs.txt"
# titre_echant <- "/work/adanguy/these/pipeline/020820/intervals/intervals.txt"        
# titre_csre <-  "/work/adanguy/these/pipeline/020820/sources/csre_genetic_map.txt"   
# titre_sortie_temp <- "/work/adanguy/these/pipeline/020820/gini/pour_gini.txt"             
# titre_graphe <-"/work/adanguy/these/pipeline/020820/graphs/gini_populations.png"    
# titre_sortie<- "/work/adanguy/these/pipeline/020820/gini/gini_populations.txt" 

cat("\n\n Input 1 : Summary PHASE output \n\n")
res <- fread(titre_resume)
head(res)
res <- res %>% filter(w_center==T) %>% na.omit() %>%
  mutate(d=lambda_rho_med*(lpop*1e6)) %>%
  dplyr::select(population, chr, region, posSNPlpop, posSNPrpop, d)
  

cat("\n\n Input 2 : intervals \n\n")
e <- fread(titre_echant) 
head(e)
e <- e %>%
  dplyr::select(chr, intervalle, posSNPlint, posSNPrint) %>%
  unique() 

cat("\n\n Input 3 : CsRe genetic map \n\n")
csre <- fread(titre_csre)
head(csre)
res <- csre %>% 
  dplyr::select(population, chr, region, posSNPlpop, posSNPrpop, dbay) %>%
  rename(d=dbay) %>%
  rbind(res)

chr <- "3A"
r <- "R1"

for (chr in unique(res$chr)){
  
  print(chr)
  
  h1 <- res %>% filter(chr==!!chr)
  
  
  for (r in unique(res$region)){
    
    
    h2 <- h1 %>% filter(region==!!r)
    
    
    minmax <- h2 %>%
      summarise(minimum=min(posSNPlpop), maximum=max(posSNPrpop))
    
    
    
    
    echant2 <- e %>% filter(chr==!!chr & posSNPrint <= !!minmax$maximum & posSNPlint >= !!minmax$minimum )
    
    
    
    
    htemp <- h2 %>%
      full_join(echant2, by=c("chr")) %>%
      filter(posSNPlint >= posSNPlpop & posSNPrint <= posSNPrpop) %>%
      dplyr::select(population, chr, region, intervalle, d)
    
    
    if (chr == unique(res$chr)[1] & r == unique(res$region)[1]){
      
      write.table(htemp, titre_sortie_temp, col.names = T, row.names = F, dec=".", sep="\t", quote=F, append = F)
      
      
      
    } else {
      
      write.table(htemp, titre_sortie_temp, col.names = F, row.names = F, dec=".", sep="\t", quote=F, append = T)
      
      
    }
    #h2 <- rbind(h2, htemp)
    
    
    
  }
  
  
  
  
  
  
  
  
  
}

chr="1A"
r="R1"

 
    
# Par chrregion
    
refEA <- fread(titre_sortie_temp) %>%
      pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
      na.omit() %>%
  group_by(chr, region) %>%
      arrange(desc(EA)) %>%
      mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
      mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
      dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
      pivot_longer(cols=c(WA,WE,EE, CsRe), names_to = "pop")  %>%
      mutate(popref="EA") %>%
      rename(ref=EA)
    
  

refEE <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  group_by(chr, region) %>%
  arrange(desc(EE)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WA,WE,EA, CsRe), names_to = "pop")  %>%
  mutate(popref="EE") %>%
  rename(ref=EE)



refWE <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  group_by(chr, region) %>%
  arrange(desc(WE)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WA,EE,EA,CsRe), names_to = "pop") %>%
  mutate(popref="WE") %>%
  rename(ref=WE)



refCsRe <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  group_by(chr, region) %>%
  arrange(desc(CsRe)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WA,EE,EA,WE), names_to = "pop") %>%
  mutate(popref="CsRe") %>%
  rename(ref=CsRe)


refWA <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  group_by(chr, region) %>%
  arrange(desc(WA)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WE,EE,EA, CsRe), names_to = "pop")  %>%
  mutate(popref="WA")%>%
  rename(ref=WA)


fusion <- rbind(refEA, refEE, refWE, refWA, refCsRe) %>%
  ungroup() %>%
  mutate(pop=factor(pop, levels = c("WE","EE","WA","EA","CsRe"))) %>%
  mutate(popref=factor(popref, levels = c("WE","EE","WA","EA","CsRe"))) %>%
  mutate(region=factor(region, levels = c("R1","R2a","C","R2b","R3"))) %>%
  group_by(chr, region, popref, pop) %>%
  mutate(gini=1-2*auc(ref, value))  
  


sortie <- fusion %>% dplyr::select(chr, region, popref, pop, gini) %>%
  unique() %>%
  as.data.frame() %>%
  droplevels() %>%
  mutate(P1.2 = case_when(popref=="WE"|pop=="WE" ~ "WE",
                          (popref!= "WE" & pop !="WE") & (popref=="EE" | pop =="EE") ~ "EE",
                          (popref != "WE" & pop !="WE" & popref !="EE" & pop != "EE") & (popref=="WA"|pop=="WA") ~ "WA",
                          (popref != "WE" & pop !="WE" & popref !="EE" & pop != "EE" & popref !="WA" & pop != "WA") & (popref=="EA"|pop=="EA")~ "EA",
                          (popref != "WE" & pop !="WE" & popref !="EE" & pop != "EE" & popref !="WA" & pop != "WA" & popref !="EA"|pop!="EA") & (popref=="CsRe"|pop=="CsRe")~ "CsRe")) %>%
  mutate(P2.2 = ifelse(popref==P1.2, as.character(pop), as.character(popref))) %>%
  dplyr::select(-one_of("popref","pop")) %>%
  rename(P1=P1.2, P2=P2.2) %>%
  group_by(P1,P2, chr, region) %>%
  summarise(gini=mean(gini)) %>%
  mutate(chrregion=paste0(chr,region)) %>%
  dplyr::select(chrregion, P1,P2, gini)


  

# genome



refEA <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  arrange(desc(EA)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WA,WE,EE, CsRe), names_to = "pop")  %>%
  mutate(popref="EA") %>%
  rename(ref=EA)



refEE <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  arrange(desc(EE)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WA,WE,EA, CsRe), names_to = "pop")  %>%
  mutate(popref="EE") %>%
  rename(ref=EE)



refWE <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  arrange(desc(WE)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WA,EE,EA,CsRe), names_to = "pop") %>%
  mutate(popref="WE") %>%
  rename(ref=WE)



refCsRe <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  arrange(desc(CsRe)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WA,EE,EA,WE), names_to = "pop") %>%
  mutate(popref="CsRe") %>%
  rename(ref=CsRe)
head(refCsRe)

refWA <- fread(titre_sortie_temp) %>%
  pivot_wider(id_cols = c("chr","region","intervalle"), values_from = "d", names_from = "population") %>%
  na.omit() %>%
  arrange(desc(WA)) %>%
  mutate(EA=cumsum(EA), WA=cumsum(WA), WE=cumsum(WE), EE=cumsum(EE), CsRe=cumsum(CsRe)) %>%
  mutate(EA=EA/max(EA), WA=WA/max(WA), WE=WE/max(WE), EE=EE/max(EE), CsRe=CsRe/max(CsRe)) %>%
  dplyr::select(intervalle, EA, WA, WE, EE, CsRe) %>%
  pivot_longer(cols=c(WE,EE,EA, CsRe), names_to = "pop")  %>%
  mutate(popref="WA")%>%
  rename(ref=WA)


fusion <- rbind(refEA, refEE, refWE, refWA, refCsRe) %>%
  ungroup() %>%
  mutate(pop=factor(pop, levels = c("WE","EE","WA","EA","CsRe"))) %>%
  mutate(popref=factor(popref, levels = c("WE","EE","WA","EA","CsRe"))) %>%
  group_by(popref, pop) %>%
  mutate(gini=1-2*auc(ref, value))  

data <- fusion %>% ungroup() %>%
  group_by(popref, pop) %>%
  arrange(value) %>% filter(ref >= 0.8) %>%
  slice(1) %>%
  ungroup() 

mod <- lm(value~gini, data=data)
tab <- data.frame(gini=seq(0,1,0.1), ref=rep(0.8, times=11), value=predict(mod, data.frame(gini=seq(0,1,0.1))))
tab

texte <- fusion %>% dplyr::select(popref, pop, gini) %>%
  mutate(gini=round(gini, digits=2)) %>%
  unique() %>%
  as.data.frame() %>%
  droplevels()
  

texte
unique(fusion$popref)

  
g <- fusion %>%  
  droplevels() %>%
    ggplot(aes(x = ref, y=value)) + geom_line() + 
  geom_abline(slope=1, intercept=0) +
  geom_ribbon(aes(ymin=value,ymax=ref), fill="blue") +
  facet_grid(popref~pop) +
  geom_text(
    size    = 5,
    data    = texte,
    mapping = aes(x = 0, y = 0.5, label = gini),
    hjust   = 0,
    vjust   = 0,
    col="blue") +
  theme_light() +
  scale_x_continuous(limits = c(0,1)) + 
  scale_y_continuous(limits=(c(0,1))) +
  xlab("% of total historical genetic distance") + ylab("% of total historical genetic distance")

cat("\n\n Graph 1 : Gini coefficients \n\n")
ggsave(titre_graphe, g)


cat("\n\n Output 1 : Gini coefficients per chrregion \n\n")
head(sortie)
write.table(sortie, titre_sortie, col.names = T, row.names = F, dec=".", sep="\t", quote=F)




sessionInfo()
