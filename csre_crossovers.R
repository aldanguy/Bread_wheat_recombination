
Sys.time()
rm(list = ls())
graphics.off()
set.seed(1)
variables <- commandArgs(trailingOnly=TRUE)


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(zoo))
suppressPackageStartupMessages(library(ggpubr))


cat("\n\ncsre_crossovers.R\n\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")



titre_genotyping <-variables[1]
titre_scaffolds <- variables[2]
titre_csre_genetic_map <-variables[3]
titre_csre_co_location <- variables[4]
iteration <- as.numeric(variables[5])
titre_position_co <- variables[6]

#  titre_genotyping <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/matrice_genotypage_csre.txt"
#  titre_scaffolds <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/scaffolds.txt"
#  titre_csre_genetic_map <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_genetic_map.txt"       
#  titre_csre_co_location <- "/work/adanguy/these/pipeline/030420/csre/co_aleatoires_csre.txt"        
#  iteration <- "1000"                                                                   
# titre_scaffolds <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/scaffolds.txt"



# titre_csre_genetic_map <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/csre_genetic_map.txt"
# titre_csre_co_location<- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/co_aleatoires.txt"
# iteration <- 10
# titre_genotyping <-"/home/adanguydesd/Documents/These_Alice/pipeline/160320/matrice_genotypage_csre.txt"
# titre_scaffolds <-"/home/adanguydesd/Documents/These_Alice/pipeline/160320/scaffolds.txt"



# Coordinates of scaffolds bourndaries
scaffolds <- fread(titre_scaffolds, header=T, dec=".", sep="\t", data.table = F)


# Genotyping matrix
# Individuals are named as X...
#######################################################################################################
matrice <-  fread(titre_genotyping, header=T, sep="\t", dec=".") 


##### Step 1 : imputation of missing data

# Missing data on monotone segments within chromosome are imputated : AA--AA--AA -> AAAAAAAA

matrice[matrice=="-"] <- NA # Imputation function catch NA values


# Forward
matriceFromLastT <- matrice %>%
  group_by(chr) %>%
  arrange(chr, posSNP) %>%
  apply(2, function(x) na.locf0(x, fromLast=T)) %>% # Case AA--BB -> AAAABB
  as.data.frame()


# Backward
matriceFromLastF <- matrice %>%
  group_by(chr) %>%
  arrange(chr, posSNP) %>%
  apply(2, function(x) na.locf0(x, fromLast=F)) %>% # Case AA--BB -> AABBBB
  as.data.frame()

# Supress letters when backward and forward desagree 
# Case AA--BB for example
matriceFromLastT[is.na(matriceFromLastT)] <- matriceFromLastF[is.na(matriceFromLastT)]
matriceFromLastT <- as.matrix(matriceFromLastT)
matriceFromLastT[matriceFromLastT != matriceFromLastF] <- "-"

matrice1 <- matriceFromLastT %>% as.data.frame() %>%
  arrange(chr, posSNP) %>%
  ungroup()
  

# No missing data remaining, except when there is a parental allele switch on a segment  : AA--BB



##### Step 2 : convert genotyping matrix in strings of parental alleles



matrice2 <- matrice1 %>%
  group_by(chr) %>%
  rename(SNPlpop=SNP, posSNPlpop=posSNP) %>%
  mutate(SNPrpop=lead(SNPlpop), posSNPrpop=lead(posSNPlpop)) %>%
  dplyr::select(population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, everything()) %>%
  mutate_at(vars(starts_with("X")), funs(paste0(lag(.), ., lead(.), lead(., n=2)))) %>% # for each individual and each SNP, paste four adjacent surrounding SNP as a string : example AABB
  na.omit() %>%
  mutate_at(vars(starts_with("X")), funs(str_remove_all(., pattern = "NA"))) # remove NA in strings

# do not worry with warnings


matrice2[1:10,1:9] %>% as.data.frame()


# 
# middle <- matrice2 %>% 
#   group_by(chr) %>%
#   slice(2:(n()-1)) %>% # remove chromosome extremities in a first place
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = c("AABB|BBAA"), replacement = "1"))) %>% # SNP whith CO with no missing data
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'B')), 4))), collapse = "|"), replacement= "0"))) %>% # all other locations 
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'B','-')), 4))), collapse = "|"), replacement= "NA")))

#middle[1:10,1:9]


# top <- matrice2 %>% 
#   group_by(chr) %>%
#   slice(1) %>%  # String patterns are not the same in chromosome extremities. left extremities here
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = c("ABB|BAA"), replacement= "1"))) %>% # SNP whith CO with no missing data
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'B')), 3))), collapse = "|"), replacement= "0"))) %>%
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'B','-')), 3))), collapse = "|"), replacement= "NA"))) 
# 
# 
# 
# 
# down <- matrice2 %>% 
#   group_by(chr) %>%
#   slice(n()) %>%  # right extremities here
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = c("AAB|BBA"), replacement= "1"))) %>%
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'B')), 3))), collapse = "|"), replacement= "0"))) %>%
#   mutate_at(vars(starts_with("X")), funs(str_replace_all(., pattern = paste0(do.call(paste0,expand.grid(rep(list(c('A', 'B','-')), 3))), collapse = "|"), replacement= "NA"))) 

##### Step 3 : CO counting






csre <- matrice1 %>% dplyr::select(population, chr, region, SNP, posSNP) %>%
  mutate(posSNP=as.numeric(as.character(posSNP))) %>%
  group_by(chr) %>%
  mutate(ligne=row_number()) %>%
  ungroup() 


# 
# marey3B <- data.frame(chr="3B", posl=seq(0,1e9,1e6), posr=seq(1e6+1,1e9+1e6+1,1e6))
# head(marey3B)

c=6
chr=as.factor("1B")


co <- data.frame()
autre <- data.frame()
chromosomes <- as.factor(as.character(unique(matrice1$chr)))
nbcotot <- 0

for (chr in chromosomes){ # loop on chromosomes
  
  print(chr)
  
  
  
  
  for (c in grep("X", names(matrice1))){ # loop on individuals
    

  ind = names(matrice1)[c] # name of individual
  
  y=as.character(matrice1[which(matrice1$chr %in% chr),c]) # string of a chromosome of an individual
  y2 <- paste0(y, collapse = "")
  nbNAy <- unlist(strsplit(y2, split="A"))
  nbNAy <- unlist(strsplit(nbNAy, split="B"))
  nbNAy <- nbNAy[which(nbNAy !="")]
  nbNAy <- unique(sort(nchar(nbNAy))) 
  nbNAy <- c(0, nbNAy) # range of missing data between parental switches
  i=0

 
    
    while(i < length(nbNAy)){ # loop on increasing missing data between parental switches
      
      i=i+1
      n=nbNAy[i]
      
      #print(n)
      

      
      
      
      
      
      # catch of intervals with crossovers
      lignes_min <- c(which(rollapply(y, 4+n, identical, c("A","A",rep("-", times=n), "B","B"))), which(rollapply(y, 4+n, identical, rev(c("A","A",rep("-", times=n), "B","B"))))) +1
      
      # Take into account chromosome extremities (less confidence because strings would look like ABB or BAA, so possibility of genotyping error)
      lignes_left <- c(which(rollapply(y, 3+n, identical, c("A",rep("-", times=n), "B","B"))), which(rollapply(y, 3+n, identical, rev(c("A","A",rep("-", times=n), "B")))))
      lignes_left <- min(lignes_left)
      lignes_right <- c(which(rollapply(y, 3+n, identical, c("A","A",rep("-", times=n), "B"))), which(rollapply(y, 3+n, identical, rev(c("A",rep("-", times=n), "B","B")))))  +1
      lignes_right <- max(lignes_right)
      
      
      
      # If a CO happend in chromosomic left extremity
      if (lignes_left ==1){
        
        lignes_min <- c(lignes_min, lignes_left)
        
      }
      
      # If a CO happend in chromosomic right extremity
      
      if (lignes_right + n + 1 ==length(y)){
        
        lignes_min <- c(lignes_min, lignes_right)
        
      }
      
      
      
      
      nbco <- length(lignes_min)
      if (nbco >=1){
        
        lignes_min <- sort(lignes_min)
        
        
        lignes_max <- lignes_min + n +1
        

        
        IDco <- paste0("C0", (nbcotot+1):(nbcotot+nbco) )
        
        
        temp <- data.frame(lignes_min=lignes_min, lignes_max=lignes_max) %>%
          mutate(chr=chr, IDco=IDco, ind=ind, nbNA=n) %>%
          left_join(csre, by="chr") %>%
          filter(ligne>=lignes_min & ligne <= lignes_max) %>%
          group_by(IDco) %>%
          rename(SNPlpop=SNP, posSNPlpop=posSNP) %>%
          mutate(SNPrpop=lead(SNPlpop), posSNPrpop=lead(posSNPlpop))%>%
          na.omit() %>%
          mutate(L=max(posSNPrpop)-min(posSNPlpop)) %>%
          rowwise() %>%
          mutate(l=posSNPrpop-posSNPlpop) %>%
          mutate(L=L/1e6) %>%
          mutate(l=l/1e6) %>%
          mutate(proba=l/L) %>% # probability to have the CO is proportional to the physical size of the interval
          ungroup() %>%
          dplyr::select(chr, IDco, ind, nbNA, SNPlpop, SNPrpop, posSNPrpop, posSNPlpop, proba, L ) %>%
          as.data.frame()
        
        # lignes <- which(matrice2$SNPlpop=="AX-89434634" & matrice2$SNPrpop=="AX-89696369")
        # matrice2[seq(lignes-1,lignes+1,1),ind]
        
       # temp2 <- temp %>% group_by(chr, IDco, ind, nbNA) %>%
        #summarise(posl=min(posSNPlpop), posr=max(posSNPrpop)) 
         
        
        
        co <- rbind(co, temp)
        # autre <- rbind(autre, temp2)
        
        nbcotot <- nbcotot + nbco
        
      }
      
      
    }
    
  }
  
  
}  
  

























cat("\n\n total nbCO \n\n")
print(nbcotot)
cat("\n\n")

##### Step 4 : CO location


nind <- length(grep("X", names(matrice1)))

co2 <- co %>%
  group_by(IDco) %>%
  mutate(probasup=cumsum(proba)) %>%
  mutate(probainf=lag(probasup)+1e-9) %>%
  ungroup() %>%
  mutate(probainf=ifelse(is.na(probainf), 0, probainf)) %>%
  na.omit() %>%
  dplyr::select(chr, IDco, L, posSNPlpop, posSNPrpop, SNPlpop, SNPrpop, proba, probainf, probasup) # range of probability to have a co at each interval














co3 <- data.frame()

for (i in 1:iteration){ # sample location of CO 
  
  print(i)

attribution <- co2 %>%
  group_by(IDco) %>%
  mutate(alea=runif(1)) %>%
  ungroup() %>%
  mutate(y=ifelse(alea >= probainf & alea <= probasup, 1, 0)) %>%
  group_by(chr, SNPlpop, SNPrpop) %>%
  summarise(y=sum(y)) %>%
  mutate(iteration=i) %>%
  dplyr::select(chr, SNPlpop, SNPrpop, y, iteration) %>%
  as.data.frame()


co3 <- rbind(co3, attribution)
  
 

}


csre3 <- matrice2 %>%
  dplyr::select(population, chr, region, SNPlpop,SNPrpop, posSNPlpop, posSNPrpop) %>%
  mutate(posSNPlpop=as.numeric(as.character(posSNPlpop))) %>%
  mutate(posSNPrpop=as.numeric(as.character(posSNPrpop))) %>%
  group_by(chr) %>%
  bind_rows(replicate(iteration-1, ., simplify = FALSE)) %>%
  group_by(SNPlpop, SNPrpop) %>%
  mutate(iteration=seq(1,iteration,1)) %>%
  full_join(co3, by=c("chr","SNPlpop","SNPrpop","iteration" )) %>%
  mutate(y=ifelse(is.na(y),0,y)) %>%
  ungroup() %>%
  mutate(M=nind) %>%
  dplyr::select(SNPlpop, SNPrpop, iteration, y, M)
  

head(csre3)


write.table(csre3, titre_csre_co_location, sep="\t", col.names = T, row.names = F, dec=".", quote=F)




  

# csre2 <- matrice2 %>% dplyr::select(population, chr, region, SNPlpop,SNPrpop, posSNPlpop, posSNPrpop) %>%
#   mutate(posSNPlpop=as.numeric(as.character(posSNPlpop))) %>%
#   mutate(posSNPrpop=as.numeric(as.character(posSNPrpop))) %>%
#   full_join(scaffolds, by="chr") %>%
#   filter(!(posSNPlpop >= posrscaf | posSNPrpop <= poslscaf)) %>%
#   group_by(population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop) %>%
#   summarise(nscaf=n()) %>%  # catch intervals whose SNP are located on two different scaffolds
#   unique() %>%
#   mutate(different_scaffold=ifelse(nscaf==1, F, T)) %>%
#   left_join(co, by=c("chr", "SNPlpop","SNPrpop")) %>%
#   group_by(population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop, different_scaffold) %>%
#   summarise(y=sum(proba, na.rm=T)) %>%
#   mutate(M=nind) %>%
#   mutate(lpop=posSNPrpop-posSNPlpop) %>%
#   mutate(lpop=lpop/1e6) %>%
#   mutate(methode="RILs") %>%
#   arrange(chr, posSNPlpop) %>%
#   ungroup() %>%
#   dplyr::select(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop,different_scaffold, lpop, y, M)
# 


csre2 <- matrice2 %>% dplyr::select(population, chr, region, SNPlpop,SNPrpop, posSNPlpop, posSNPrpop) %>%
  mutate(posSNPlpop=as.numeric(as.character(posSNPlpop))) %>%
  mutate(posSNPrpop=as.numeric(as.character(posSNPrpop))) %>%
  full_join(scaffolds, by="chr") %>%
  filter(!(posSNPlpop >= posrscaf | posSNPrpop <= poslscaf)) %>%
  group_by(population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop) %>%
  summarise(nscaf=n()) %>%  # catch intervals whose SNP are located on two different scaffolds
  unique() %>%
  mutate(different_scaffold=ifelse(nscaf==1, F, T)) %>%
  mutate(lpop=posSNPrpop-posSNPlpop) %>%
  mutate(lpop=lpop/1e6) %>%
  mutate(methode="RILs") %>%
  arrange(chr, posSNPlpop) %>%
  ungroup() %>%
  dplyr::select(methode, population, chr, region, SNPlpop, SNPrpop, posSNPlpop, posSNPrpop,different_scaffold, lpop)





head(csre2)


write.table(csre2, titre_csre_genetic_map, sep="\t", col.names = T, row.names = F, dec=".", quote=F)




write.table(co, titre_position_co, sep="\t", col.names = T, row.names = F, dec=".", quote=F)




















sessionInfo()

