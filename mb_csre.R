
Sys.time()
rm(list = ls())
graphics.off()
set.seed(1)
variables <- commandArgs(trailingOnly=TRUE)


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))


cat("\n\nmb_csre.R\n\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")



titre_co_position <-variables[1]
titre_csre <- variables[2]
titre_chr_tab <- variables[3]
iteration <- as.numeric(variables[4])
titre_fenetres <- variables[5]
titre_sortie <- variables[6]
# 
# titre_co_position <- "/work/adanguy/these/pipeline/030420/csre/position_co.txt"  
# titre_csre <-  "/work/adanguy/these/pipeline/030420/sources/csre_genetic_map.txt" 
# titre_chr_tab <- "/work/adanguy/these/pipeline/amont/Decoupage_chr_ble.tab"         
# iteration <- as.numeric("10")                                                         
# titre_fenetres <- "/work/adanguy/these/pipeline/030420/temp/carte_chrregion.txt"     
# titre_sortie <- "/work/adanguy/these/pipeline/030420/csre/carte_chrregion_csre.txt"



cat("\n\n Input 1 : Random assignment of CO \n\n")
co <- fread(titre_co_position)
head(co)

cat("\n\n Input 2 : CsRe genetic map \n\n")
csre <- fread(titre_csre)
head(csre)

cat("\n\n Input 3 : windows \n\n")
carte1M <- fread(titre_fenetres)
head(carte1M)



haldane <- function(prop_rec_RILs){
  
  prop_rec_gam <- prop_rec_RILs/(2-2*prop_rec_RILs)
  
  return(prop_rec_gam)
}


# minimum <- 1
# maximum <- 1e9


regions  <- read.table(titre_chr_tab, header=T, dec=".", sep="\t", skip=1) %>%
  dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
  mutate(chr = str_remove(Chromosome,"chr")) %>%
  dplyr::select(-one_of("Chromosome"))


shape <- csre %>% dplyr::select(region, shape, rate) %>%
  unique()
  

minmax <- csre%>%
  group_by(chr) %>% summarise(min=min(posSNPlpop), max=max(posSNPrpop))


nind <- unique(csre$M)


carte1M <- carte1M %>% full_join(minmax, by="chr") %>%
  filter(!(posl >= max | posr <= min)) %>%
  dplyr::select(chr, posl, posr) %>%
  arrange(chr, posl)

co <- co %>%
  group_by(IDco) %>%
  mutate(probasup=cumsum(proba)) %>%
  mutate(probainf=lag(probasup)+1e-9) %>%
  ungroup() %>%
  mutate(probainf=ifelse(is.na(probainf), 0, probainf)) %>%
  na.omit() %>%
  dplyr::select(chr, IDco, L, posSNPlpop, posSNPrpop) # range of probability to have a co at each interval




coM <- co %>%
  group_by(IDco) %>%
  full_join(carte1M, by="chr") %>%
  filter(! (posSNPlpop >= posr | posSNPrpop <= posl )) %>%
  na.omit() %>%
  group_by(IDco,posSNPlpop, posSNPrpop, posl, posr) %>%
  mutate(L1=min(posSNPrpop, posr) - max(posSNPlpop, posl) ) %>% 
  mutate(L=L*1e6) %>% # L = window of the CO
  mutate(proba=L1/L) %>% # % of the CO window overlapping the interval of 1Mb
  group_by(IDco) %>%
  mutate(probasup=cumsum(proba)) %>%
  mutate(probainf=lag(probasup)+1e-9) %>%
  mutate(probainf=ifelse(is.na(probainf), 0, probainf)) %>%
  ungroup() %>%
  mutate(L=L/1e6) %>%
  dplyr::select(chr, IDco, L, posl, posr, proba, probainf, probasup) 




co3 <- data.frame()

for (i in 1:iteration){ # sample location of CO 
  
  print(i)
  
  attribution <- coM %>%
    group_by(IDco) %>%
    mutate(alea=runif(1)) %>%
    ungroup() %>%
    mutate(y=ifelse(alea >= probainf & alea <= probasup, 1, 0)) %>%
    group_by(chr, posl, posr) %>%
    summarise(y=sum(y)) %>%
    mutate(iteration=i) %>%
    dplyr::select(chr, posl, posr, y, iteration) %>%
    as.data.frame()
  
  
  co3 <- rbind(co3, attribution)
  
  
  
}



co4 <- carte1M %>% 
  group_by(chr) %>%
  bind_rows(replicate(iteration-1, ., simplify = FALSE)) %>%
  group_by(chr, posl, posr) %>%
  mutate(iteration=seq(1,iteration,1)) %>%
  full_join(co3, by=c("chr","posl","posr","iteration" )) %>%
  mutate(y=ifelse(is.na(y),0,y)) %>%
  ungroup() %>%
  mutate(M=nind)%>%
  full_join(regions, by="chr") %>%
  mutate(region = ifelse(posl <=  R1.R2a*1e6, "R1",
                         ifelse(posl >= R1.R2a*1e6 +1 & posl <= R2a.C*1e6, "R2a",
                                ifelse(posl >= R2a.C*1e6 +1 & posl <= C.R2b*1e6, "C", 
                                       ifelse(posl >= C.R2b*1e6 +1 & posl <= R2b.R3*1e6, "R2b",
                                              ifelse(posl >= R2b.R3*1e6 +1, "R3", NA)))))) %>%
  full_join(shape, by="region") %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  arrange(chr, posl, iteration, region) %>%
  group_by(chr, region, posl, posr, iteration,y, shape, rate) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(l=(posr-posl)/1e6) %>%
  mutate(Cbay=(y+shape)/((M*l) +rate)) %>%
  mutate(sdCbay=(y+shape)/(((M*l) +rate)^2)) %>%
  mutate(Dbay=Cbay*l) %>%
  mutate(Rbay=Dbay) %>%
  mutate(rbay=haldane(Rbay)) %>%
  mutate(dbay=rbay) %>%
  mutate(cbay=dbay/l) %>%
  mutate(cbay=100*cbay) %>%
  group_by(chr, region, posl, posr, shape, rate, l) %>%
  summarise(cbay=mean(cbay), sdCbay=mean(log10(sdCbay))) %>% ## ou sdCbay=mean(sdCbay)
  mutate(population="CsRe") %>%
  arrange(chr, posl) %>%
  ungroup() %>%
  dplyr::select(population, chr, region, posl, posr, cbay, sdCbay)




cat("\n\n Output 1 : CsRe genetic map at 1 Mb \n\n")
head(co4)  
write.table(co4, titre_sortie, col.names = T, row.names = F, dec=".", sep="\t", quote=F)

# 
# g <- co4 %>% filter(chr=="3B") %>%
#   ggplot(aes(x=posl, y=cbay, col=region)) + geom_line() +
#   ylab("Bayesian rec rates (cM/Mb)") +
#   xlab("physical position (pb)") +
#   ggtitle("3B") +
#   theme_light() +
#   theme(plot.title = element_text(hjust = 0.5))
# g
# 
# cat("\n\n Graph 1 : CsRe genetic map at 1 Mb for chr 3B \n\n")
# ggsave(titre_graphe, g, width = 7.29, height = 4.5, units = c("in"))
# 

sessionInfo()