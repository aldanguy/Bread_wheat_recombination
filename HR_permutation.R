




Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bedr))




cat("\n HR_permutation.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

titre_resume <-variables[1]
titre_HR <- variables[2]
titre_chr_tab <- variables[3]
seed <- as.numeric(variables[4])
titre_permutation <- variables[5]


# titre_resume <- "/work/adanguy/these/pipeline/020820/PHASE/PHASE_summary_outputs.txt"    
# titre_HR <-  "/work/adanguy/these/pipeline/020820/HR/HR.txt"                          
# titre_chr_tab <-  "/work/adanguy/these/pipeline/amont/Decoupage_chr_ble.tab"               
# seed <- as.numeric("621")                                                                
# titre_permutation <- "/work/adanguy/these/pipeline/020820/HR/permutations/permutation_621.txt"


set.seed(seed)
nb_essai_max <- 1e4

  

cat("\n\n Input 1 : HR\n\n")
h <- fread(titre_HR)
head(h)
h <- h  %>%
  arrange(population, chr, posl) %>%
  filter(kept==T)  %>%
  unique()

cat("\n\n Input 2 : summary output of PHASE \n\n")
res <- fread(titre_resume)
head(res)
res <-res  %>%
  dplyr::select(chr, region, posSNPlpop, posSNPrpop)



cat("\n\n Input 3: limits of region \n\n")
region <- read.table(titre_chr_tab, header=T, dec=".", sep="\t", skip=1)
head(region)
regionsr  <- region %>%
  dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
  mutate(chr = str_remove(Chromosome,"chr")) %>%
  dplyr::select(-one_of("Chromosome")) %>%
  rename(R1=R1.R2a, R2a=R2a.C, C=C.R2b, R2b=R2b.R3) %>%
  mutate(R2a=R2a*1e6, C=C*1e6, R2b=R2b*1e6, R1=R1*1e6) %>%
  pivot_longer(cols=c("R1","R2a","C","R2b"), names_to = "region", values_to = "posr")


regionsl  <- region %>%
  dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
  mutate(chr = str_remove(Chromosome,"chr")) %>%
  dplyr::select(-one_of("Chromosome")) %>%
  rename(R2a=R1.R2a, C=R2a.C, R2b=C.R2b, R3=R2b.R3) %>%
  mutate(R2a=R2a*1e6 +1, C=C*1e6 +1, R2b=R2b*1e6 +1, R3=R3*1e6 +1) %>%
  pivot_longer(cols=c("R2a","C","R2b","R3"), names_to = "region", values_to = "posl")



regionsr <- res %>% group_by(chr, region) %>%
  summarise(posr=max(posSNPrpop)) %>%
  filter(region=="R3") %>%
  ungroup() %>%
  rbind(.,regionsr)


regionsl <- res %>% group_by(chr, region) %>%
  summarise(posl=min(posSNPlpop)) %>%
  filter(region=="R1") %>%
  ungroup() %>%
  rbind(.,regionsl)

region <- regionsr %>% full_join(regionsl, by=c("chr","region")) %>%
  arrange(chr, posl) %>%
  dplyr::select(chr, region, posl, posr) %>%
  as.data.frame() 



rm(res, regionsl, regionsr)


chr1 <- unique(h$chr)[1]
r1 <- unique(h$region)[1]
p1=unique(h$population[1])



chr="3A"
p="EA"
r="R1"
permutationf <- data.frame()



for (chr in unique(h$chr)){
  
  print(chr)
  
  
  h1 <- h %>% filter(chr==!!chr)
  
  
  
  
  
  
  for (r in unique(h$region)){
    
    

    
    
    h2 <- h1 %>% filter(region==!!r)
    
    
    minmax <- region %>% filter(chr==!!chr & region==!!r) 
    
    
    pos_possible <- seq(minmax$posl, minmax$posr,1)
    
    
    for (p in unique(h$population)){
      
      

      
      
      h3 <- h2 %>% filter(population==!!p)
      
      
      
      
      
      
      
      
      
      
      
      
      
      if (nrow(h3)>=1){
        
        reussite=F
        k=0
        
        
        
       # while(reussite==F & k<=nb_essai_max){
          while(reussite==F){
            
          
          if(k==nb_essai_max){
            
            cat("\n\n probleme ; error when attributing new position to HR \n\n")
            print(p)
            print(chr)
            print(r)
            cat("\n\n")
            stop()
            
            
          }
          
          
          posl <- sample(pos_possible,nrow(h3), replace=F)
          posr <- posl + h3$taille 
          pos <- bedr.snm.region(paste0("chr1:",posl,"-",posr),  verbose = F, check.chr = F, check.valid = F, check.zero.based=F)
        

          
          
          
          if (length(pos)==nrow(h3)) {
            
            
            permutation <- data.frame(population=p,
                                      chr=chr,
                                      region=r,
                                      posl=posl,
                                      posr=posr,
                                      taille=h3$taille,
                                      IDintHR=h3$IDintHR,
                                      kept=T, 
                                      lambda_med=NA,
                                      d_genet_histo=NA,
                                      iteration=seed)
            
            

            
            
            ################# quelle forme pour que ça ressemble à HR de début coloc
            

            
            reussite <- T
            
            
            if (chr == chr1 & p ==  p1 & r == r1){
              cat("\n\n Output 1 : Permutation of HRs\n\n")
              print(head(permutation))
              write.table(permutation, titre_permutation, col.names = T, row.names = F, dec=".", sep="\t", quote=F, append = F)
              
              
            } else {
              
              write.table(permutation, titre_permutation, col.names = F, row.names = F, dec=".", sep="\t", quote=F, append = T)
              
              
              
            }
            
          }
          
          k=k+1
          
          
        }
        
        
        
        
        
        
        
        
      }
      
      
    }
    
    
  }
  
}


sessionInfo()
