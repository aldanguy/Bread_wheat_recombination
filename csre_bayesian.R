

Sys.time()
rm(list = ls())
graphics.off()
set.seed(1)
variables <- commandArgs(trailingOnly=TRUE)


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))






cat("\n\ncsre_bayesian.R\n\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_csre_intervals <- variables[1]
titre_csre_co_location <- variables[2]
dossier_graphs <- variables[3]
titre_ecdf <- variables[4]
titre_mks <- variables[5]

 # titre_csre_intervals <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_intervals.txt"
 # titre_csre_co_location <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/anciens/co_aleatoires_csre.txt"
 # dossier_graphs <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/graphes/"
 # titre_mks <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_mks.txt"
 # titre_csre <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_genetic_map.txt"


# titre_csre_intervals <- "/work/adanguy/these/pipeline/030420/sources/csre_genetic_map.txt"
# titre_csre_co_location <- "/work/adanguy/these/pipeline/030420/csre/co_aleatoires_csre.txt" 
# dossier_graphs <- "/work/adanguy/these/pipeline/030420/graphes/"                    
# titre_mks <- "/work/adanguy/these/pipeline/030420/temp/mqs/csre_mqs.txt"  


# Formula from Haldane Waddington 1931
# Convert the proportion of recombined RILs to the proportion of recombined gametes

haldane <- function(prop_rec_RILs){
  
  prop_rec_gam <- prop_rec_RILs/(2-2*prop_rec_RILs)
  
  return(prop_rec_gam)
}





# Empirical estimates of RILs recombination rates per interval
cat("\n\n Input 1 : Intervals in CsRe population and frequency of CO per interval \n\n")
# y = number of CO interpolated if missing data
# M = number of RILs
csre <- fread(titre_csre_intervals, header=T, dec=".", sep="\t", data.table = F) 
head(csre)




cat("\n\n Input 2 : Allocation of CO to CsRe intervals according probabilities \n\n")
# Previously, we identify CO : AA--BB or BB--AA. Parental switch were separated by a various number of missing data (-).
# Each CO was located in a window L (A--...--B). L had a physical length.
# Each interval overlapped by L had its own physical size l
# The ratio l/L indicate the probabilty to receive the CO in this interval (hypothesis : rec rate constant across L)
# We sample CO according l/L for each intervals in 1 000 independant iteration
# y = sum of CO in one interval in one iteration
co <- fread(titre_csre_co_location, header=T, dec=".", sep="\t", data.table = F) 
head(co)
# Big file (> 6Gb)
# use filter(iteration <= 5)



# Previously, we identify CO : AA--BB or BB--AA. Parental switch were separated by a various number of missing data (-).
# Each CO was located in a window L (A--...--B). L had a physical length.
# Each interval overlapped by L had its own physical size l
# The ratio l/L indicate the probabilty to receive the CO in this interval (hypothesis : rec rate constant across L)
# y = sum of probability in one interval

csre2 <- csre %>%
  dplyr::select(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold, lpop) %>%
  mutate(region=factor(region, levels = c("R1","R2a","C","R2b","R3"))) %>%
  full_join(co, by=c("SNPlpop","SNPrpop"), suffix=c("",".x")) %>%
  na.omit() %>%
  mutate(Remp=y/M) %>% # proportion of recombined RILs
  mutate(Demp=Remp) %>% # Genetic distance of RILs population. Approximation D= R (no Haldane or Kosambi mapping function)
  mutate(Cemp=Demp/lpop) %>% # RILs recombination rates
  mutate(Cemp=100*Cemp) %>% # Conversion in cM/Mb
  mutate(remp=haldane(Remp)) %>% # proportion of recombined gametes
  mutate(demp=remp) %>% # meiotic genetic distance (expected to be > 50 cM and < 600 cM = between 1 and 3 CO/meiosis/chr)
  mutate(cemp=demp/lpop) %>% # meiotic rec rates
  mutate(cemp=100*cemp) %>%  # Conversion in cM/Mb
  group_by(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold, lpop,M) %>%
  summarise(ymean=mean(y),
            ymax=max(y),
            Remp=mean(Remp),
            Demp=mean(Demp),
            Cemp=mean(Cemp),
            remp=mean(remp),
            demp=mean(demp),
            cemp=mean(cemp)) %>%
  ungroup() %>%
  dplyr::select(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold, lpop,M,ymean,ymax,Remp,Demp,Cemp,remp,demp,cemp)
  
cat("\n\n cumulated genetic distance per chr (M) \n\n")
csre2 %>% group_by(chr) %>% summarise(sum(demp)) %>% as.data.frame()

nbinter0 <- csre2 %>% filter(Cemp ==0) %>% nrow(.)
nbinter <- nrow(csre2)

cat("\n\n prop inter rec=0 \n\n")
print(round(nbinter0/nbinter, digits=2))
cat("\n\n")

png(titre_ecdf)
plot(x=log10(csre2$cemp), y=ecdf(csre2$cemp)(csre2$cemp), ylim=c(0,1), xlab="emprical recombination rates (cM/Mb)", ylab="ECDF")
dev.off()

# Position of CO randomly placed among possible locations over 1 000 iterations
# So nb_rows = 1 000 iterations * nb of intervals in CsRe map
# co <- fread(titre_csre_co_location, header=T, dec=".", sep="\t", data.table = F)


##### Step 1 : prior distribution of recomibation rates in RILs


prior_gamma <- data.frame()
i=0
b="R1"
g <- list()
for (b in c("R1","R3","R2a","R2b","C")){
  
  
  print(b)
  
  

  # empirical estimates of RILs recombination rates in region b
  C <- csre2$Cemp[which(csre$region==b)]
  
  # number of datapoints to fit a Gamma distribution
  print(length(C))
  
  
  
  # nNll recombination rates are replace by 1e-5 M/Mb, which seems low enough
  #C[C==0] <- 0.001
  C[C==0] <- min(C[C>0])
  print(min(C))
  
  # Conversion of RILs recombination rates from cM/Mb to M/Mb
  # It is importante that Gamma parameters units are in M instead of cM
  # Because R = y/M is in M unit
  C <- C*0.01
  #hist(log10(C), main=b, xlab="recombination rates observed in RILs (M/Mb)")

  
  
  # package MASS_7.3-51.4
  ajuste <- fitdistr(C,"gamma", method = c("Nelder-Mead"))
  shape <-  ajuste$estimate[1]
  rate <- ajuste$estimate[2]
  
  # A tab to save Gamma parameters value
  i=i+1
  prior_gamma[i,"region"] <- b
  prior_gamma[i,"shape"] <- shape
  prior_gamma[i,"rate"] <- rate
 
  
  # Simulations of Gamma, to draw the distribution
  simu <- rgamma(1e6, shape = shape, rate = rate)
 
   # A list of plots, one per region
   g[[i]] <- ggplot(data=data.frame(C), aes(x=C)) +
     geom_density()  +
     geom_density(data=data.frame(simu), aes(x=simu), col="red") +
     scale_x_continuous(limits = c(1e-5,0.01))+
     scale_y_continuous(limits = c(0, 3100)) +
     geom_hline(yintercept=0, colour="white", size=1) +
     ggtitle(b) +
     theme(plot.title = element_text(hjust = 0.5)) +
     xlab("recombination rates observed in RILs (M/Mb)")  +
     annotate("text", x=Inf, y = Inf, vjust=1, hjust=1 ,label =  paste0("Gamma (shape = ", round(shape,digits=2), ", rate = ",round(rate,digits=2),")"), col="red")
   
   
   rm(simu)
   rm(C)

   
   

}

cat("\n\n prior paramaters \n\n")
prior_gamma

# Combined graph for prior
p0 <- do.call(grid.arrange,g)
titre_graphe_p0 <- paste0(dossier_graphs,"prior.png")
cat("\n\n Graph 1 : prior and empirical distributions of RILs recombination rates per region \n\n")
ggsave(titre_graphe_p0,p0)
rm(p0)
rm(g)


csre$shape <- prior_gamma$shape[match(csre$region, prior_gamma$region)]
csre$rate <- prior_gamma$rate[match(csre$region, prior_gamma$region)]


##### Step 2 : compute the mean of the posterior distribution of CsRe RILs recombination rates


co <- csre %>%
  dplyr::select(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold, lpop, shape, rate) %>%
  mutate(region=factor(region, levels = c("R1","R2a","C","R2b","R3"))) %>%
  full_join(co, by=c("SNPlpop","SNPrpop"), suffix=c("",".x")) %>%
  na.omit() %>%
  mutate(Cbay=(y + shape)/(M*lpop + rate)) %>% # formula of conjugate bayesian model. RILs recombination rates, in M/Mb
  mutate(Dbay=Cbay*lpop) %>% # genetic distance of RILs population
  mutate(Rbay=Dbay) %>% # theorical proportion of recombined RILs according to bayesian model. Approximation R=D
  mutate(rbay=haldane(Rbay)) %>% # theorical proportion of recombined gametes
  mutate(dbay=rbay) %>% # Approximation r=d
  mutate(cbay=dbay/lpop) %>% # bayesian estimates of meiotic recombination rates
  mutate(Cbay=100*Cbay) %>% # convert in cM/Mb
  mutate(cbay=100*cbay) %>% # convert in cM/Mb
  group_by(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold, lpop,M, shape, rate) %>%
  summarise(Rbay=mean(Rbay),
            Dbay=mean(Dbay),
            Cbay=mean(Cbay),
            rbay=mean(rbay),
            dbay=mean(dbay),
            cbay=mean(cbay)) %>%
  ungroup() %>%
  dplyr::select(SNPlpop, SNPrpop,Rbay,Dbay,Cbay,rbay,dbay,cbay, shape, rate)

co <- co %>% full_join(csre2, by=c("SNPlpop", "SNPrpop")) %>%
  mutate(ygraphe= round(ifelse(ymean >=5, 5, ymean)))%>% # y_max is the maximal number of CO per interval observed accross the 1000 iteration
  mutate(ygraphe=as.character(ygraphe)) %>%
  mutate(ygraphe=ifelse(ygraphe=="5","5+",ygraphe)) %>%
  dplyr::select(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold, lpop,M,ymean, ymax,ygraphe, Remp,Rbay,Demp,Dbay,Cemp,Cbay,remp,rbay,demp,dbay,cemp, cbay, shape, rate)


# 
# co <- csre2 %>%
#   mutate(Cbay=(y + shape)/(M*lpop + rate)) %>% # formula of conjugate bayesian model. RILs recombination rates, in M/Mb
#   mutate(Dbay=Cbay*lpop) %>% # genetic distance of RILs population
#   mutate(Rbay=Dbay) %>% # theorical proportion of recombined RILs according to bayesian model. Approximation R=D
#   mutate(rbay=haldane(Rbay)) %>% # theorical proportion of recombined gametes
#   mutate(dbay=rbay) %>% # Approximation r=d
#   mutate(cbay=dbay/lpop) %>% # bayesian estimates of meiotic recombination rates
#   mutate(Cbay=100*Cbay) %>% # convert in cM/Mb
#   mutate(cbay=100*cbay) %>% # convert in cM/Mb
#   group_by(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold, lpop, y, M, Cemp, cemp, shape, rate )%>%
#   dplyr::summarise(Cbay=mean(Cbay), y_mean=mean(y.x), y_max = max(y.x), cbay=mean(cbay)) %>% # summarise over the 1000 sample of the number of CO per intervals
#   rowwise() %>%
#   arrange(chr, posSNPlpop) %>%
#   mutate(demp=cemp*lpop) %>%
#   mutate(dbay=cbay*lpop) %>%
#   mutate(y_max = min(y_max, 5)) %>%
#   mutate(y_max=as.character(y_max)) %>%
#   mutate(y_max= ifelse(y_max=="5","5+",y_max))%>% # y_max is the maximal number of CO per interval observed accross the 1000 iteration
#   dplyr::select(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop,different_scaffold, lpop, M, Cemp, cemp, demp, Cbay, cbay, dbay, shape, rate, y, y_mean,y_max  )

  


# head(co) %>% as.data.frame()

cat("\n\n effect of bayesian model on the cumulated genetic distance per chromosome \n\n")
co %>% mutate(ratio=cbay/cemp) %>% group_by(chr) %>%
  summarise(dbay=sum(dbay), demp=sum(demp)) %>% ungroup() %>%
  mutate(ratio=dbay/demp) %>% mutate(ratiomean=round(median(ratio), digits=2), ratiosd=round(IQR(ratio), digits=2)) %>% as.data.frame()


# bottom left triangle = intervals with very low probability(=0) to observe a CO in empirical method, and no CO was observe during 1000 iterations
# higher bayesian estimates = lower l and over-estimated y (and a lower y too, but related to the fact tha a lower l decrease the probability to observe a co)
# plateau at cbay ~ 2.4 = when lpop = 1e-6 (can not be lower)
# co %>% filter(region=="R1" & cemp <2.5 & cbay <2.5 & cemp >0 & cbay >2 & y_max==1) %>% as.data.frame()
# co <- co %>% mutate(y2=ifelse(y>=0.1,y,0)) %>%
#   mutate(Remp2=y2/406) %>%
#   mutate(remp2=haldane(Remp2)) %>%
#   mutate(demp2=remp2) %>%
#   mutate(cemp2=demp2/lpop) %>%
#   mutate(cemp2=cemp2*100) %>%
#   mutate(ratio=y2/y_mean) 
# co %>% filter(ratio >=3 & ratio < Inf) %>% head()
# 
#  co %>%
#   ggplot(aes(x=cemp, y=cbay, colour=ratio)) +
#   geom_point(size=1) +
#   scale_x_continuous(limits=c(0,10)) +
#   scale_y_continuous(limits=c(0,10)) +
#   geom_abline(slope=1, intercept=0) +
#   xlab("empirical recombination rates (cM/Mb)") +
#   ylab("bayesian recombination rates (cM/Mb)") +
#   theme_light()
# value beyond the triangle
#co %>% filter(region=="R1" & cemp <2.5 & cbay <5 & cemp >0 & cbay >2.5) %>% as.data.frame()
# y_mean is overestimated y (with true estimate of y, these points should be lower )


##### Step 4 : graphs


p1 <- co %>%
  ggplot(aes(x=cemp, y=cbay, colour=ygraphe)) +
  geom_point(size=1) +
  scale_x_continuous(limits=c(0,10)) +
  scale_y_continuous(limits=c(0,10)) +
  scale_color_brewer(palette = "Paired") +
  geom_abline(slope=1, intercept=0) +
  guides(color=guide_legend(title="# CO")) +
  xlab("empirical recombination rates (cM/Mb)") +
  ylab("Bayesian recombination rates (cM/Mb)") +
  theme_light()
titre_graphe_p1 <- paste0(dossier_graphs,"bayesian.png")
cat("\n\n Graph 2 : Empirical and bayesian estimates of meiotic recombinnation rates of CsRe (colors function of maximum nb of observed CO over the 1000 iteration ) \n\n")
ggsave(titre_graphe_p1, p1, width = 7.29, height = 4.5, units = c("in"))



p2 <- co %>% 
  filter(region=="R1") %>%
  filter((ymean>0.07 & ymean < 0.8) | (ymean > 0.01 & ymean < 0.02)) %>%
  ggplot(aes(x=cemp, y=cbay, colour=log10(ymean))) +
  geom_point(size=1) +
  scale_x_continuous(limits=c(0,2.5)) +
  scale_y_continuous(limits=c(0,2.5)) +
  geom_abline(slope=1, intercept=0) +
  xlab("empirical recombination rates (cM/Mb)") +
  ylab("Bayesian recombination rates (cM/Mb)") +
  theme_light() +
  facet_wrap(.~region)

titre_graphe_p2 <- paste0(dossier_graphs,"bayesian_y.png")
cat("\n\n Graph 2 : Empirical and bayesian estimates of meiotic recombination rates of CsRe (colors function of CO frenquency) \n\n")
ggsave(titre_graphe_p2, p2, width = 7.29, height = 4.5, units = c("in"))



p3 <- co %>% 
  ggplot(aes(x=cemp, y=cbay, colour=log10(ymean))) +
  geom_point(size=0.005) +
  scale_x_continuous(limits=c(0,10)) +
  scale_y_continuous(limits=c(0,10)) +
  geom_abline(slope=1, intercept=0) +
  xlab("empirical recombination rates (cM/Mb)") +
  ylab("Bayesian recombination rates (cM/Mb)") +
  theme_light() +
  facet_wrap(.~region)


titre_graphe_p3 <- paste0(dossier_graphs,"bayesian_region.png")
cat("\n\n Graph 3 : Empirical and bayesian estimates of meiotic recombinnation rates of CsRe (colors function of physical size of the interval ) \n\n")
ggsave(titre_graphe_p3, p3, width = 7.29, height = 4.5, units = c("in"))


cat("\n\n average cM/Mb \n\n")
co %>% group_by(chr) %>% summarise(dg=100*sum(dbay), dp=sum(lpop)) %>%
  mutate(r=dg/dp) %>%
  ungroup() %>%
  summarise(r=mean(r))


co %>% ungroup() %>%
  mutate(extreme=case_when(cemp==0~"null",
                                  cemp >=100 ~ "high")) %>%
  mutate(ntot=n()) %>%
  na.omit() %>%
  group_by(extreme, ntot) %>%
  summarise(n2=n()) %>%
  mutate(n3=n2/ntot)


##### Step 5 : polymorphic markers in CsRe


mqs <- co %>% dplyr::select(SNPlpop, SNPrpop) %>%
  gather() %>%
  dplyr::select(value) %>%
  unique() %>%
  unlist() %>%
  as.vector()




cat("\n\n Output 1 : Empirical and bayesian estimates of meiotic recombinnation rates of CsRe \n\n")
write.table(co, titre_csre_intervals, col.names = T, row.names = F, dec=".", sep="\t", quote = F)
head(co) %>% as.data.frame()
cat("\n\n Output 2 : polymorphic markers in CsRe \n\n")
write.table(mqs, titre_mks, col.names = F, row.names = F, dec=".", sep="\t", quote = F)
head(mqs)

sessionInfo()

