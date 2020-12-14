


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




cat("\n gini_genetic_distance.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_resume <-variables[1]
titre_csre <- variables[2]
titre_HR <- variables[3]
titre_graphe <- variables[4]
titre_sortie<- variables[5]


   # # titre_graphe <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/graphes/gini.png"
   # titre_resume <-"/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/PHASE_summary_outputs.txt"
   # titre_csre <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/csre_genetic_map.txt"
   # titre_HR <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/HR.txt"

cat("\n\n Input 1 : Summary PHASE output \n\n")
res <- fread(titre_resume)
head(res)
res <- res %>% filter(w_center==T) %>% na.omit() %>%
  mutate(lpop=lpop*1e6) %>%
  mutate(d=lambda_rho_med*(lpop)) %>%
  rename(l=lpop) %>%
  dplyr::select(population, chr, region, posSNPlpop, posSNPrpop, l, d)
  



cat("\n\n Input 2 : CsRe genetic map \n\n")
csre <- fread(titre_csre)
head(csre)
res <- csre %>% 
  dplyr::select(population, chr, region, posSNPlpop, posSNPrpop, dbay, lpop) %>%
  mutate(lpop=lpop*1e6) %>%
  rename(d=dbay, l=lpop) %>%
  rbind(res)

cat("\n\n Input 3 : HR \n\n")
hr <- fread(titre_HR)
head(hr)
hr <- hr %>% filter(kept==T)

  

# genome




res2 <- res %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA","CsRe"))) %>%
  group_by(population) %>%
  arrange(population, desc(d)) %>%
  mutate(d2=cumsum(d)/sum(d)) %>%
  mutate(l2=cumsum(l)/sum(l))


texte <- res2 %>%
  group_by(population) %>%
  summarise(gini=round(1-2*unique(auc(d2,l2)), digits=2)) %>%
  ungroup() %>%
  mutate(mean_gini=mean(gini))
texte

d_pop <- res2 %>%  group_by(population) %>%
  filter(population != "CsRe") %>%
  filter(d2 >= 0.8) %>% slice(1) %>%
  dplyr::select(population, d2, l2) %>%
  rename(prop_genet_d=d2, prop_d_phy=l2) %>%
  mutate(prop_d_phy=round(prop_d_phy, digits = 2)) %>%
  ungroup() %>%
  mutate(mean_prop_d_phy=round(mean(prop_d_phy), digits=2), sd_propr_d_phy=round(sd(prop_d_phy), digits=2))

d_csre <- res2 %>%  group_by(population) %>%
  filter(population == "CsRe") %>%
  filter(d2 >= 0.8) %>% slice(1) %>%
  dplyr::select(population, d2, l2) %>%
  rename(prop_genet_d=d2, prop_d_phy=l2) %>%
  mutate(prop_d_phy=round(prop_d_phy, digits = 2)) %>%
  ungroup()

d_pop
d_csre


texte <- d_pop %>% dplyr::select(population, prop_genet_d, prop_d_phy) %>%
  rbind(., d_csre) %>%
  inner_join(texte, by="population") %>%
  mutate(prop_genet_d=round(prop_genet_d,2))

color <- as.character(res2$population)
color[which(color=="WE")] <- "#D7191C"
color[which(color=="EE")] <- "#FDAE61"
color[which(color=="WA")] <- "#ABDDA4"
color[which(color=="EA")] <- "#2B83BA"
color[which(color=="CsRe")] <- "yellow"


graphe <- res2 %>%
  ggplot(aes(y=l2,x=d2)) + geom_line() +
  theme_light() +
  geom_abline(intercept=0, slope=1) +
  geom_ribbon(aes(ymin=d2,ymax=l2), fill=color) +
  facet_wrap(population~., ncol=2) +
  ylab("cumulated physical distance (%)") +
  xlab("cumulated genetic distance (%)") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        strip.text.x = element_text(size = 14, color="black"),
        axis.text.y = element_text(size=14),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.title.x = element_text(size=16),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        strip.background = element_rect(color="black", fill="white"),
        panel.grid.minor = element_blank())+
  geom_segment(data = texte, aes(x=0.8, xend=0.8, y=0, yend= prop_d_phy), col="black", linetype="dashed") +
  geom_segment(data = texte, aes(x=0, xend=0.8, y=prop_d_phy, yend= prop_d_phy), col="black", linetype="dashed") +
  geom_text(
    size    = 5,
    data    = texte,
    mapping = aes(x = 0, y = 0.9, label = paste0("Gini ",gini)),
    hjust   = 0,
    vjust   = 0,
    col="blue") +
  geom_text(
    size    = 5,
    data    = texte,
    mapping = aes(x = 0, y = 0.5, label = prop_d_phy),
    hjust   = 0,
    vjust   = 0,
    col="black") +
  geom_text(
    size    = 5,
    data    = texte,
    mapping = aes(x = 0.85, y = 0, label = prop_genet_d),
    hjust   = 0,
    vjust   = 0,
    col="black")+  scale_x_continuous(limits = c(0,1), breaks=c(0,0.25,0.5,0.8,1)) +  scale_y_continuous(limits = c(0,1))
graphe

cat("\n\n Graph 1 : Proportion of genetic distance in physical distance \n")
ggsave(titre_graphe, graphe)



d <- res2 %>%  group_by(population) %>%
  filter(population != "CsRe") %>%
  summarise(dtot=sum(d), ltot=sum(l)) 

cat("\n\n Statistics on HR \n")
hr %>% group_by(population) %>%
  summarise(d=sum(d_genet_histo), l=sum(taille)) %>%
  inner_join(d, by="population") %>%
  mutate(ratio_d=round(d/dtot, digits=2), ratio_l=round(l/ltot, digits=3), sd(ratio_l))

 
sessionInfo()
