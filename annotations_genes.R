

Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))




cat("\n annotations_genes.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

titre_annotation <- variables[1]
titre_HR <- variables[2]
titre_chr_tab <- variables[3]
titre_resume <- variables[4]
titre_sortie <- variables[5]
titre_graphe <- variables[6]

# titre_resume <-"/home/adanguydesd/Documents/These_Alice/pipeline/160320/sorties_resumees_PHASE.txt"
# titre_annotation <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/annotation.txt"
# titre_HR <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/hrec.txt"
# titre_chr_tab <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Decoupage_chr_ble.tab"
# titre_sortie <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/annotation_hotspots.txt"
# titre_graphe <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/graphes/annotation_hotspots.png"

cat("\n\n Input 1 : annotation\n\n")
a <- fread(titre_annotation, header=F)
head(a)
a <- a %>% rename(chr=V1, feature=V2, posl=V3,posr=V4, ID=V5) %>%
  filter(feature %in% c("gene", "exon", "five_prime_UTR", "three_prime_UTR")) %>%
  mutate(feature=factor(feature, levels=c("gene","five_prime_UTR","exon","three_prime_UTR")))%>%
  arrange(feature, chr, posl, posr) %>%
  filter(! row_number() %in% which(duplicated(.[c("chr","posl","posr")]))) %>%
  group_by(ID) %>%
  mutate(posmin=min(posl), posmax=max(posr)) %>%
  filter(! (feature !="gene" & posl==posmin & posr ==posmax)) %>%
  ungroup() %>%
  dplyr::select(chr, feature, posl, posr, ID) %>%
  arrange(ID, posl, desc(posr)) 


cat("\n\n Input 2 Hotspots\n\n")
h <- fread(titre_HR)
head(h)
h <- h  %>%
  arrange(population, chr, posl) %>%
  filter(kept==T) %>%
  mutate(hr=T) %>%
  dplyr::select(population, chr, posl, posr, IDintHR) %>%
  unique()

a %>% filter(feature=="gene") %>%
  mutate(taille=posr-posl) %>%
  summarise(mean_size_gene=mean(taille), sd_size_gene=sd(taille))


cat("\n\n Input 3 Summary PHASE output \n\n")
res <- fread(titre_resume)
head(res)
h <- res  %>%
  na.omit() %>%
  filter(w_center==T) %>%
  dplyr::select(population, chr, region, posSNPlpop, posSNPrpop) %>%
  full_join(h, by=c("chr","population","posSNPlpop"="posl")) %>%
  filter(population=="WE") %>%
  mutate(hr=ifelse(is.na(IDintHR),F,T)) %>%
  mutate(posl=posSNPlpop) %>%
  mutate(posr=posSNPrpop) %>%
  dplyr::select(population, chr, region, posl, posr, IDintHR, hr) %>%
  arrange(chr, posl) 
  
  
  


cat("\n\n Input 3 regions \n\n")
regions  <- read.table(titre_chr_tab, header=T, dec=".", sep="\t", skip=1)
head(regions)
regions <- regions %>%
  dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
  mutate(chr = str_remove(Chromosome,"chr")) %>%
  dplyr::select(-one_of("Chromosome"))



a <- a %>% full_join(regions, by='chr') %>%
  mutate(region = ifelse(posl <=  R1.R2a*1e6, "R1",
                         ifelse(posl >= R1.R2a*1e6 +1 & posl <= R2a.C*1e6, "R2a",
                                ifelse(posl >= R2a.C*1e6 +1 & posl <= C.R2b*1e6, "C", 
                                       ifelse(posl >= C.R2b*1e6 +1 & posl <= R2b.R3*1e6, "R2b",
                                              ifelse(posl >= R2b.R3*1e6 +1, "R3", NA)))))) %>%
  dplyr::select(chr, region, feature, posl, posr, ID)

chr <- "3A"
r <- "R1"
p <- "WE"
f=unique(a$feature)[1]


for (chr in unique(h$chr)) {
  
  print(chr)
  
  
  h1 <- h %>% filter(chr==!!chr)
  
  for (p in unique(h$population)){
    print(p)
    
    
    h2 <- h1 %>% filter(population==!!p)
    
    
    for (f in unique(a$feature)){
      
      
      for (r in unique(h$region)){
        
        h3 <- h2 %>% filter(region==!!r)
        
        
        
        
        # h2 <- h %>% filter(population==!!p & chr==!!chr & region==!!r) %>%
        #   full_join(a2, by="chr") %>%
        #   filter(!(posr <= posSNPlpop | posl >= posSNPrpop)) %>%
        #   rename(region=region.x) %>%
        #   group_by(population, chr, region, feature,pc) %>%
        #   summarise(n=n()) %>%
        #   pivot_wider(values_from = "n", names_from = "pc") %>%
        #   mutate(ratio=`TRUE`/(`TRUE`+`FALSE`))
        # 
        
        
        
        minmax <- h3 %>%
          summarise(minimum=min(posl), maximum=max(posr))
        
        a2 <- a %>% filter(chr==!!chr & feature==!!f & posl <= !!minmax$maximum & posr >= !!minmax$minimum )
        
        
        nTassociated <- h3 %>% 
          filter(hr==T) %>%
          full_join(a2, by="chr") %>%
          filter(!(posr.y <= posl.x | posl.y >= posr.x)) %>%
          dplyr::select(posl.x) %>%
          unique() %>%
          nrow()
        
        
        
        
        nFassociated <- h3 %>%
          filter(hr==F) %>%
          full_join(a2, by="chr") %>%
          filter(!(posr.y <= posl.x | posl.y >= posr.x)) %>%
          dplyr::select(posl.x) %>%
          unique() %>%
          nrow()
        
        
        
        
        tab <- data.frame(population=p,
                          chr=chr,
                          region=r,
                          feature=f,
                          nTassociated=nTassociated,
                          nFassociated=nFassociated,
                          nTtot=nrow(h %>% filter(population==!!p & chr==!!chr & region==!!r & hr==T)),
                          nFtot=nrow(h %>% filter(population==!!p & chr==!!chr & region==!!r & hr==F)))
        
        
        if (p==unique(h$population)[1] & chr==unique(h$chr)[1] & r==unique(h$region)[1] & f==unique(a$feature)[1]){
          
          
          write.table(tab, titre_sortie, row.names = F, dec=".", sep="\t", quote=F, col.names = T, append = F)
          
        } else {
          
          write.table(tab, titre_sortie, row.names = F, dec=".", sep="\t", quote=F, col.names = F, append = T)
        }
        
      }
      
    }
    
  }
  
  
}

tab <- fread(titre_sortie)

cat("\n\n Association of hotspots and features")
head(tab)
write.table(tab, titre_sortie, row.names = F, dec=".", sep="\t", quote=F, col.names = T, append = F)




tab <- tab %>% 
  group_by(region, feature) %>%
  summarise(nTassociated=sum(nTassociated), nTtot=sum(nTtot), nFassociated=sum(nFassociated), nFtot=sum(nFtot) ) %>%
  mutate(Tassociated = nTassociated/nTtot, Fassociated=nFassociated/nFtot, Tnonassociated=1-Tassociated, Fnonassociated=1-Fassociated) %>%
  dplyr::select(region, feature, Tassociated, Fassociated) %>%
  rename(HR=Tassociated, non_HR=Fassociated) %>%
  pivot_longer(-c("region","feature"), names_to = "HR", values_to = "per" ) %>%
  ungroup() %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) 

tab <- tab %>% filter(feature=="gene") %>%
  mutate(per=1-per, feature="intergenic") %>%
  rbind(., tab)  %>%
  mutate(feature=factor(feature, levels=c("intergenic", "gene","five_prime_UTR","exon","three_prime_UTR")))
write.table(tab, titre_sortie, row.names = F, dec=".", sep="\t", quote=F, col.names = T, append = F)


tab %>% group_by(feature, HR) %>%
  summarise(mean_prop=round(mean(per, na.rm=T), digits=2), sd_prop=round(sd(per, na.rm=T), digits=2))

g <- tab %>%
  mutate(HR=ifelse(HR=="HR","TRUE","FALSE")) %>%
  mutate(HR=factor(HR, levels=c("TRUE","FALSE"))) %>%
  ggplot(aes(x=feature, y=per, fill=HR)) + geom_bar(stat = "identity", color="black", position=position_dodge())+
  theme_light() +
  ylab("Percentage of intervals overlapping gene feature") +
  xlab("") +
  facet_grid(region~.)





ggsave(titre_graphe,g)


sessionInfo()