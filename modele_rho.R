



### Compute correlations of historical recombination profiles (lambda or rho) of 4 populations of landraces per chr*region with a heterocedastic mixed model
### Outputs are : variance-covariance matrix (*aleatoires*) ; residuals (*ajuste*) ; fixed effets (*fixes*)
## Inputs are : one estimate of lambda (or rho) per population and per interval. 


# module purge
# module load system/R-3.4.3_bis

Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(asreml))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gtools))




cat("\n modele_rho.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_data <- variables[1]
titre_aleatoires_chrregion <- variables[2]
titre_fixes_chrregion <- variables[3]
titre_ajuste_chrregion <- variables[4]

# titre_data <- "/work/adanguy/these/pipeline/030420/fusion/echant_de_novo_intervals_balanced.txt"
# titre_aleatoires_chrregion <- "/work/adanguy/these/pipeline/030420/mod/chrregion_lambda_aleatoires.txt"         
# titre_fixes_chrregion <-"/work/adanguy/these/pipeline/030420/mod/chrregion_lambda_fixes.txt"              
# titre_ajuste_chrregion <-  "/work/adanguy/these/pipeline/030420/mod/chrregion_lambda_ajuste.txt"             


# titre_data <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/echant_de_novo_intervals_balanced.txt"



########## Fonctions

makeSymm1 <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

makeSymm2 <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

LRT <- function(L0,L1, d0, d1) { 
  LRT = abs((L1-L0))
  value=qchisq(0.95,abs(d0-d1))
  if (LRT >= value){ aleatoires <- TRUE} else {aleatoires <- FALSE}
  return(aleatoires)}


diff_AIC <- function(L0,L1, d0, d1) {
  
  AIC0 <- -2 * L0 + 2*d0
  AIC1 <- -2 * L1 + 2*d1
  
  dAIC = AIC0 - AIC1 > 0
  
  # expected to be TRUE if AIC0 > AIC1  
  #http://faculty.washington.edu/skalski/classes/QERM597/papers_xtra/Burnham%20and%20Anderson.pdf
  return(dAIC)
  
  
}



diff_BIC <- function(L0,L1, d0, d1, nedf0, nedf1) {
  
  BIC0 <- -2 * L0 + d0 * log(nedf0)
  BIC1 <- -2 * L1 + d1 *log(nedf1)
  
  dBIC = BIC0 - BIC1 > 0
  #Kass, Robert E.; Raftery, Adrian E. (1995), "Bayes Factors", Journal of the American Statistical Association, 90 (430): 773–795, doi:10.2307/2291091, ISSN 0162-1459, JSTOR 2291091.
  return(dBIC)
  
  
}

#################### Donnees


cat("\n\n Input 1 : Recombination rates estimates for 5 populations \n\n")
datar <- fread(titre_data, header=T, sep="\t", dec=".", data.table = F)
head(datar)



datar <- datar %>%
   dplyr::select(population, chr, region, intervalle, SNPlint, SNPrint, posSNPlint,posSNPrint, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop,lint,lpop, different_scaffold, ID, lambda, recombinaison ) %>%
   arrange(population, chr,posSNPlint) %>%
   mutate(different_scaffold = as.factor(different_scaffold)) %>%
   mutate(intervalle = factor(intervalle, levels=unique(datar$intervalle))) %>%
   mutate(chr=as.factor(chr)) %>%
   mutate(region = factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
   mutate(population = factor(population, levels=c("CsRe","EA","EE","WA","WE"))) %>%
   mutate(chrregion=as.factor(paste0(chr,region)))%>%
   arrange(population, chr, posSNPlint) %>%
   dplyr::select(population, chr, region, chrregion, intervalle, SNPlint, SNPrint, posSNPlint,posSNPrint, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, lint, lpop, different_scaffold, ID, lambda, recombinaison )

datar <- datar %>%
  mutate(chrregion=factor(chrregion, levels=as.character(unique(datar$chrregion)))) %>%
  arrange(chrregion, posSNPlint, population)


datar <- datar %>% filter(population !="CsRe")  %>%
  rename(lambda_rho=recombinaison) %>%
  mutate(y=log10(lambda_rho))  %>%
  droplevels()


summary(datar)
head(datar)

pop=as.character(rev(unique(datar$population)))
npop=length(pop)
nb_cases=npop + npop*(npop-1)/2
liste_chrregion <- as.character(unique(datar$chrregion))


cat("\n\n genomewide correlation\n\n")
datar %>% dplyr::select(population, intervalle, y) %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
  pivot_wider(id_cols=intervalle, values_from = y, names_from = population) %>%
  dplyr::select(-intervalle) %>%
  cor() %>%
  round(digits=2)

chrregion="7DR3"
for (chrregion in liste_chrregion){
  
  cat("\n\nchrregion")
  print(chrregion)
  cat("\n\n")
  
  
  datat <- datar %>% filter(chrregion==!!chrregion) %>% droplevels() %>% arrange(population, intervalle)
  

  
  nb_intervalles <- length(unique(datat$intervalle))
  
  # y = log10(rho/rec_csre)
  modas1 <- asreml(fixed = y ~ population,
                   rcov=~corh(population):intervalle,
                   data=datat,
                   maxiter=1000)
  
  
  
  modas2 <- asreml(fixed = y ~ population,
                   rcov=~corgh(population):intervalle,
                   data=datat,
                   maxiter=1000)
  
  

  
  
  
  if (modas2$last.message=="LogLikelihood Converged"){
    
    
    L1 <- summary(modas1)$loglik
    df1 <- length(summary(modas1)$varcomp[,1])
    nedf1=summary(modas1)$nedf
    
    L2 <- summary(modas2)$loglik
    df2 <- length(summary(modas2)$varcomp[,1])
    nedf2=summary(modas2)$nedf
    
    
    LRTtest <- LRT(L1,L2,df1,df2)
    AICtest <- diff_AIC(L1, L2, df1, df2)
    BICtest <- diff_BIC(L1, L2, df1, df2, nedf1, nedf2)
    

    
    
    fixes <- modas2$coeff$fixed %>%
      as.data.frame(col.names="estimateur", fix.empty.names=T) %>%
      rownames_to_column("param") %>%
      rename_at( 2, ~"estimateur" ) %>%
      mutate(param=gsub("population_","", param)) %>%
      mutate(param=gsub("\\(Intercept\\)","mean", param)) %>%
      mutate(param=as.factor(param)) %>%
      mutate(chrregion=chrregion) %>%
      dplyr::select(chrregion, param, estimateur)
    
    
    
    param_alea = matrix(0, npop,npop)
    param_alea [lower.tri(param_alea, diag=F)]=rev(summary(modas2)$varcomp[2:((npop*(npop-1)/2) +1),1])
    param_alea <- makeSymm2(param_alea)
    diag(param_alea) <- rev(summary(modas2)$varcomp[((npop*(npop-1)/2) +2):((npop*(npop-1)/2) +1 +npop),1])
    param_alea <- as.data.frame(param_alea)
    colnames(param_alea) = row.names(param_alea) <- pop
    
    
    noms <- combinations(n = npop, r = 2, v = pop, repeats.allowed = TRUE)  
    aleatoires <- data.frame(chrregion, noms, estimateur=param_alea[noms]) %>%
      mutate(LRT=LRTtest, AIC=AICtest,BIC=BICtest, mod_complet_converge=TRUE, nb_intervalles=nb_intervalles)%>%
      rename(P1=X1, P2=X2) %>%
      mutate(nb_intervalles=nb_intervalles) %>%
      mutate(P1.2 = case_when(P1=="WE"|P2=="WE" ~ "WE",
                              (P1!= "WE" & P2 !="WE") & (P1=="EE" | P2 =="EE") ~ "EE",
                              (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE") & (P1=="WA"|P2=="WA") ~ "WA",
                              (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA") & (P2=="EA"|P2=="EA")~ "EA",
                              (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA" & P1 !="EA"& P2 !="EA" ) & (P2=="CsRe"|P2=="CsRe")~ "CsRe")) %>%
      mutate(P2.2 = ifelse(P1==P1.2, as.character(P2), as.character(P1))) %>%
      dplyr::select(-one_of("P1","P2")) %>%
      rename(P1=P1.2, P2=P2.2)  %>%
      dplyr::select(chrregion, P1, P2, estimateur, mod_complet_converge, LRT, AIC, BIC, nb_intervalles)
    
    
    
    

    
     datat2 <- datat %>%
       mutate(mean=fixes$estimateur[which(fixes$param=="mean")]) %>%
       mutate(effet_pop=fixes$estimateur[match(datat$population, fixes$param)]) %>%
       mutate(effet_pop=mean+effet_pop) %>%
       mutate(res=modas2$res) %>%
       as.data.frame()
    
    
    
  } else {
    
    

    fixes <- modas1$coeff$fixed %>%
      as.data.frame(col.names="estimateur", fix.empty.names=T) %>%
      rownames_to_column("param") %>%
      rename_at( 2, ~"estimateur" ) %>%
      mutate(param=gsub("population_","", param)) %>%
      mutate(param=gsub("\\(Intercept\\)","mean", param)) %>%
      mutate(param=as.factor(param)) %>%
      mutate(chrregion=chrregion) %>%
      dplyr::select(chrregion, param, estimateur)
    
    
    param_alea = matrix(0, npop,npop)
    param_alea [lower.tri(param_alea, diag=F)]=summary(modas1)$varcomp[2,1]
    param_alea <- makeSymm2(param_alea)
    diag(param_alea) <- summary(modas1)$varcomp[3:(npop+2),1]
    colnames(param_alea) = row.names(param_alea) <- rev(pop)
    
    
    
    
    noms <- combinations(n = npop, r = 2, v = pop, repeats.allowed = TRUE)  
    aleatoires <- data.frame(chrregion, noms, estimateur=param_alea[noms]) %>%
      mutate(LRT=NA, AIC=NA,BIC=NA, mod_complet_converge=FALSE, nb_intervalles=nb_intervalles)%>%
      rename(P1=X1, P2=X2) %>%
      mutate(nb_intervalles=nb_intervalles) %>%
      mutate(P1.2 = case_when(P1=="WE"|P2=="WE" ~ "WE",
                              (P1!= "WE" & P2 !="WE") & (P1=="EE" | P2 =="EE") ~ "EE",
                              (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE") & (P1=="WA"|P2=="WA") ~ "WA",
                              (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA") & (P2=="EA"|P2=="EA")~ "EA",
                              (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA" & P1 !="EA"& P2 !="EA" ) & (P2=="CsRe"|P2=="CsRe")~ "CsRe")) %>%
      mutate(P2.2 = ifelse(P1==P1.2, as.character(P2), as.character(P1))) %>%
      dplyr::select(-one_of("P1","P2")) %>%
      rename(P1=P1.2, P2=P2.2) %>%
      dplyr::select(chrregion, P1, P2, estimateur, mod_complet_converge, LRT, AIC, BIC, nb_intervalles)
    
    
    
    
    
    
    
    
    
    datat2 <- datat %>%
      mutate(mean=fixes$estimateur[which(fixes$param=="mean")]) %>%
      mutate(effet_pop=fixes$estimateur[match(datat$population, fixes$param)]) %>%
      mutate(effet_pop=mean+effet_pop) %>%
      mutate(res=modas1$res) %>%
      as.data.frame()
    
    
    
  }
  
  
  
  
  
  
  
  
  datat3 <- datat2 %>%   
    dplyr::select(population,
                  chr,
                  region,
                  ID,
                  intervalle,
                  SNPlint,
                  SNPlpop,
                  SNPrint,
                  SNPrpop,
                  posSNPlint,
                  posSNPlpop,
                  posSNPrint,
                  posSNPrpop,
                  lint,
                  lpop,
                  different_scaffold,
                  lambda,
                  lambda_rho,
                  res)

  
  
  
  
  
  if(chrregion==liste_chrregion[1]){
    
    

    print(head(datat2))
    
    cat("\n\n Output 1 : variance-correlation matrix parameter for the chrregion \n\n")
    write.table(aleatoires, titre_aleatoires_chrregion, col.names=T, row.names=F, dec=".", sep="\t", quote=F)
    print(head(aleatoires))
    cat("\n\n Output 2 : fixed effects estimates for the chrregion \n\n")
    write.table(fixes, titre_fixes_chrregion, col.names=T, row.names=F, dec=".", sep="\t", quote=F)
    print(head(fixes))
    cat("\n\n Output 3 : estimates of rec rates per landraces populations \n\n")
    write.table(datat3, titre_ajuste_chrregion, col.names=T, row.names=F, dec=".", sep="\t", quote=F)
    print(head(datat3))

    
    rm(fixes, datat2, datat3)
    
    
    
  } else {
    
    write.table(aleatoires, titre_aleatoires_chrregion, col.names=F, row.names=F, dec=".", sep="\t", quote=F, append=T)
    write.table(fixes, titre_fixes_chrregion, col.names=F, row.names=F, dec=".", sep="\t", quote=F, append=T)
    write.table(datat3, titre_ajuste_chrregion, col.names=F, row.names=F, dec=".", sep="\t", quote=F, append=T)

    rm(datat2, datat3)
    
    
    
    
  }
  
} 





sessionInfo()