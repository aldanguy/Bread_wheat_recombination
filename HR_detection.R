
Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bedr))




cat("\n HR_detection.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_resume <-variables[1]
titre_chr_code <- variables[2]
titre_HR <- variables[3]
titre_graphe_distri_lambda <- variables[4]
titre_graphe_choix_taille_HR <- variables[5]
titre_graphe_lambda_taille_intervalles <- variables[6]
titre_graphe_taille_HR <- variables[7]

 # titre_resume <-"/home/adanguydesd/Documents/These_Alice/pipeline/160320/sorties_resumees_PHASE.txt"
 # titre_HR <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/HR.txt"
 # titre_chr_code <-"/home/adanguydesd/Documents/These_Alice/pipeline/amont/Codes_chr.txt"

 
seuil <- 4
quantile <- 0.2
 
 
 
cat("\n\n Input 1 : code chr \n\n")
chr_code <- fread(titre_chr_code) 
chr_code
chr_code <- chr_code %>% rename(chr_code_lettre=V1, chr_code_nombre=V2)




cat("\n\n Input 2 : summary statistic of PHASE outputs \n\n")
res <- fread(titre_resume) 
head(res)
res <- res %>%
  filter(w_center==T) %>%
  na.omit() %>%
  arrange(population, chr, posSNPlpop) %>%
  inner_join(chr_code, by=c("chr"="chr_code_lettre")) %>%
  mutate(taille=posSNPrpop-posSNPlpop) %>%
  mutate(d_genet_histo=lambda_rho_med*taille) %>%
  dplyr::select(population, chr, chr_code_nombre, region, posSNPlpop, posSNPrpop, SNPlpop, SNPrpop, taille, lambda_med, d_genet_histo)






############ HRs detection


 res2 <- res %>% mutate(hr=ifelse(lambda_med>=seuil & region !="C", T, F)) %>%
   filter(hr==T)



#### Merging of adjacent HRs

 res3 <-  res2 %>%
   rename(start=posSNPlpop, stop=posSNPrpop) %>%
   dplyr::select(chr_code_nombre, start, stop, population) %>%
   rename(chr=chr_code_nombre) %>%
   arrange(population, chr, start) %>%
   mutate(chr=paste0("chr",chr)) %>%
   bedr.merge.region(., stratify.by = "population", verbose = F) %>%
   rename(population=names, posSNPlpop=start, posSNPrpop=end) %>%
   mutate(chr=as.numeric(gsub("chr","", chr))) %>%
   left_join(chr_code, by=c("chr"="chr_code_nombre")) %>%
   dplyr::select(population, chr, chr_code_lettre, posSNPlpop, posSNPrpop) %>%
   full_join(res, by=c("chr"="chr_code_nombre","population"), suffix=c("",".x")) %>%
   filter(posSNPlpop.x >= posSNPlpop & posSNPrpop.x <= posSNPrpop) %>%
   group_by(population, chr_code_lettre, region, posSNPlpop, posSNPrpop) %>%
   summarise(lambda_med=max(lambda_med), d_genet_histo=sum(d_genet_histo)) %>%
   rename(chr=chr_code_lettre) %>%
   mutate(taille=posSNPrpop - posSNPlpop) %>%
   rename(posl=posSNPlpop, posr=posSNPrpop) %>%
   ungroup() %>%
   arrange(population, chr, region, posl, posr, taille, lambda_med, d_genet_histo)


 
 
 q=0
 sortie <- data.frame()
 for (q in seq(0,1,0.01)) {
 
  res4 <- res3 %>%
    filter(taille <= quantile(taille, (1-q/2)) & taille >= quantile(taille, q/2)) %>%
    mutate(taillemax=max(taille), taillemin=min(taille), ratio=taillemax/taillemin) %>%
    dplyr::select(ratio) %>%
    mutate(quantile=q) %>%
    unique()
 
 
 # res3 <- res2 %>%
 #   ungroup() %>%
 #   mutate(taille=posSNPrpop - posSNPlpop) %>%
 #   filter(taille <= quantile(taille, (1-q/2))) %>%
 #   mutate(taillemax=max(taille), taillemin=min(taille), ratio=taillemax/taillemin) %>%
 #   ungroup() %>%
 #   dplyr::select(ratio) %>%
 #   mutate(quantile=q) %>%
 #   unique()
 
 
 sortie <- rbind(sortie, res4)
 
 
 }
 
 cat("\n\n Graph 1: histogram of lambda \n\n")
 png(titre_graphe_distri_lambda)
 hist(res$lambda_med, xlim = c(0,10), n=10000, main="", xlab=expression(lambda))
 abline(v=seuil, col="red", lty=2)
 dev.off()
 
 cat("\n\n Graph 2: variation of hotspots in a range with supressing quantile \n\n")
 png(titre_graphe_choix_taille_HR)
 plot(sortie$quantile, log10(sortie$ratio), main="", ylab="Size max/Size min (log10)", xlab="Supression of widest and smallest HR (quantile) ")
 abline(v=quantile, col="red", lty=2)
 dev.off()

 res4 <- res3 %>%
    mutate(kept=ifelse(taille <= quantile(taille, 1-(quantile/2)) & taille >= quantile(taille, (quantile/2)) ,T,F)) %>%
   group_by(population, chr, region) %>%
   mutate(IDintHR=paste(population,chr,region,seq(1,n(),1), sep="_")) %>%
   ungroup() %>%
   dplyr::select(population, chr, region, posl, posr, taille, IDintHR, kept, lambda_med, d_genet_histo)

 table(res4$kept)
 
 size <- res4 %>% filter(kept==T) %>%
   summarise(min=min(taille), max=max(taille))
 


g <- res  %>%
  mutate(HR=ifelse(lambda_med>=seuil & region !="C" & posSNPrpop-posSNPlpop <= size$max & posSNPrpop-posSNPlpop >= size$min,T,F)) %>%
  mutate(HR=factor(HR, levels=c("TRUE","FALSE"))) %>%
  ggplot(aes(x=taille, y=lambda_med, col=HR)) + geom_point(size=0.02) +
  scale_x_continuous(trans="log10") + scale_y_continuous(trans="log10") +
  xlab("Physical size (pb)") + 
  ylab(expression(lambda)) +
  geom_vline( xintercept = size$min) +
  geom_vline(xintercept = size$max) +
  theme_light()+
  geom_hline(yintercept = seuil)
  

g2 <- res %>%
  mutate(HR=ifelse(lambda_med>=seuil & region !="C" & posSNPrpop-posSNPlpop <= size$max & posSNPrpop-posSNPlpop >= size$min,T,F)) %>%
  mutate(HR=factor(HR, levels=c("TRUE","FALSE"))) %>%
  ggplot(aes(x=taille, fill=HR)) + geom_histogram(bins=100) +
  scale_x_continuous(trans="log10") +
  theme_light() +
  xlab("Physical size (pb)") +
  ylab("") +
  facet_grid(HR~., scales = "free")


cat("\n\n Graph 3: relationship between lambda and size of the intervals \n\n")
ggsave(titre_graphe_lambda_taille_intervalles, g)


cat("\n\n Graph 4: physical size of HR \n\n")
ggsave(titre_graphe_taille_HR, g2)


cat("\n\n Output 1: HR \n\n")
head(res4)
write.table(res4, titre_HR, col.names = T, row.names = F, dec=".", sep="\t", quote=F)



sessionInfo()







