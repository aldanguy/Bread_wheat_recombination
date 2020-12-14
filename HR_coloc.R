




Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(bedr))
suppressPackageStartupMessages(library(igraph))





cat("\n HR_coloc.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

type <- variables[1]
titre_HR <-variables[2]
titre_resume <- variables[3]
titre_sortie_coloc <- variables[4]
marges <- as.numeric(variables[5])
titre_sortie_cliques_HR <- variables[6]




# type="reference"
# titre_HR <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/HR.txt"
# titre_sortie_coloc <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/coloc_HR.txt"
# 
# type="100"
# titre_HR <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/100_HR.txt"
# titre_sortie_HR_100 <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/clique_100.txt"
# titre_sortie_coloc <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/coloc_100_HR.txt"
# 
# 
# type="autre"
# titre_HR <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/permutation_HR_1.txt"
# titre_sortie_coloc <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/permutation_HR_1.txt"
# 

# type= "reference"                                                            
# titre_HR <- "/work/adanguy/these/pipeline/030420/fusion/HR.txt"                  
# titre_resume <- "/work/adanguy/these/pipeline/030420/fusion/sorties_resumees_PHASE.txt"
# titre_sortie_coloc <-  "/work/adanguy/these/pipeline/030420/fusion/coloc_HR.txt"            
# marges <- 0

 # type <- "permutation"                                                    
 # titre_HR <- "/work/adanguy/these/pipeline/030420/temp/HR/permutation_1.txt"
 # titre_sortie_coloc <- "/work/adanguy/these/pipeline/030420/temp/HR/permutation_1.txt"


cat("\n\n Input 1 : HR \n\n")
h <- fread(titre_HR)
head(h)
h <- h  %>%
  filter(kept==T) %>%
  arrange(population, chr, posl) %>%
  dplyr::select(population, chr, region, posl, posr ,IDintHR, taille, lambda_med, d_genet_histo, kept) %>%
  unique()



cat("\n\n Input 2 : summary output PHASE \n\n")
res <- fread(titre_resume)
head(res)

if (is.na(marges)){
marges <- res %>% filter(w_center==T) %>% na.omit() %>%
  mutate(lpop=lpop*1e6) %>%
  group_by(region) %>% summarise(marge=median(lpop))

print(marges)

}

####### Ajdency matrix
i=1
matrice <- data.frame(matrix(0, nrow = nrow(h), ncol = nrow(h)), row.names=h$IDintHR)
colnames(matrice) <- h$IDintHR

for (i in 1: nrow(h)){
  
  marge <- marges
  
  if (!is.vector(marges)){
    
    marge <- h %>% slice(i) %>% inner_join(marges, by="region") %>% dplyr::select(marge) %>%
    unlist() %>% as.vector()
    
  }
  
  test <- h %>% slice(i) %>%
    left_join(h %>% slice(-i), by="chr") %>%
    filter(!(posl.x >= posr.y + !!marge | posr.x  <= posl.y - !!marge)) %>%
    as.data.frame()
  
  
  IDintHR.x <- h %>% slice(i) %>% dplyr::select(IDintHR) %>% unlist() %>% as.vector()
  
  IDintHR.y <- test %>% dplyr::select(IDintHR.y) %>%
    unlist() %>%
    as.vector()
  
  
  matrice[IDintHR.x, which(colnames(matrice) %in% IDintHR.y)] <- 1
  matrice[ which(colnames(matrice) %in% IDintHR.y),IDintHR.x] <- 1
  
  
  
}
#which(unlist(matrice["WE_3A_R1_20",])>0)

#### Find cliques

trouver_les_cliques <- function(matrice,h, pop_to_keep ){
  
  
  matricetemp <- matrice
  
  
  populations <- c("WE","EE","WA","EA")
  
  
  pop_to_remove <- populations[which(!populations %in% pop_to_keep)]
  
  
  i=1
  
  while(i <= length(pop_to_remove)){
    
    
    matricetemp <- matricetemp[-grep(pop_to_remove[i],row.names(matricetemp)), -grep(pop_to_remove[i],row.names(matricetemp))]
    
    i=i+1
    
  }
  
  
  
  grapĥe <- graph_from_adjacency_matrix(as.matrix(matricetemp), mode = "lower", weighted = T, diag = F, add.rownames = NA)
  test <- cliques(grapĥe, min=length(pop_to_keep), max=length(pop_to_keep))
  final <- data.frame()
  
  
  if (length(test) >=1){
    
    
    
    
    l=1
    for (l in 1:length(test)){
      
      
      pop=as.numeric(c("WE","EE","WA","EA") %in% unlist(strsplit(names(unlist(test)),split="_")))
      
      
      
      IDEA <- names(test[[l]])[grep("EA",names(test[[l]]))]
      IDEE <- names(test[[l]])[grep("EE",names(test[[l]]))]
      IDWE <-  names(test[[l]])[grep("WE",names(test[[l]]))]
      IDWA <-  names(test[[l]])[grep("WA",names(test[[l]]))]
      
      
      if (length(IDEA) ==0){
        
        IDEA <- NA
        
      } 
      
      if (length(IDWA) ==0){
        
        IDWA <- NA
        
      } 
      
      if (length(IDEE) ==0){
        
        IDEE <- NA
        
      } 
      
      if (length(IDWE) ==0){
        
        IDWE <- NA
        
      } 
      
      
      pos <- h %>% filter(IDintHR %in% !!names(test[[l]])) %>%
        summarise(posl=min(posl),posr=max(posr), chr=unique(chr))
      
      data <- data.frame(clique=l, nb = length(pop_to_keep),
                         WE=pop[1], EE=pop[2], WA=pop[3], EA=pop[4],
                         IDWE=IDWE,IDEE=IDEE,IDWA=IDWA,IDEA=IDEA,
                         posl=pos$posl,
                         posr=pos$posr,
                         chr=pos$chr,
                         coloc=paste(pop_to_keep, collapse = ""))
      
      final <- rbind(final, data)
      
    }
    
    
    final <- final %>%
      mutate(chr=as.factor(as.character(chr))) %>%
      arrange(chr, posl) %>%
      mutate(clique=seq(1,n(),1))
    
    
  }
  
  invisible(final)
}

WEEEWAEA <- trouver_les_cliques(matrice,h, c("WE","EE","WA","EA"))
WEEEWA <- trouver_les_cliques(matrice, h,c("WE","EE","WA"))
WEEEEA <- trouver_les_cliques(matrice, h,c("WE","EE","EA"))
WEWAEA <- trouver_les_cliques(matrice, h,c("WE","WA","EA"))
EEWAEA <- trouver_les_cliques(matrice, h,c("EE","WA","EA"))
WEEE <- trouver_les_cliques(matrice, h,c("WE","EE"))
WEWA <- trouver_les_cliques(matrice, h,c("WE","WA"))
WEEA <- trouver_les_cliques(matrice, h,c("WE","EA"))
EEWA <- trouver_les_cliques(matrice, h,c("EE","WA"))
EEEA <- trouver_les_cliques(matrice, h,c("EE","EA"))
WAEA <- trouver_les_cliques(matrice, h,c("WA","EA"))



final <- rbind(WEEEWAEA,
               WEEEWA,WEEEEA,WEWAEA, EEWAEA,
               WEEE,WEWA,WEEA,EEWA,EEEA, WAEA) %>%
  arrange(desc(nb), chr, posl, coloc) %>%
  mutate(n=row_number())



# test <- data.frame(clique=NA, nb=1, WE=0, EE=1,WA=1, EA=0, IDWE=NA, IDEE="EE_1A_R1_16", IDWA="WA_1A_R1_17", IDEA=NA, posr=1, posl=2, chr="1A",pop="EEWA",
#                    n=nrow(final)+1)
# final <- rbind(final, test)


i=1
for (i in 1:nrow(final)){
  
  
  
  
  
  
  to_remove <- final %>% filter(n==i) %>%
    left_join(final %>% filter(n!=i), by="chr")  %>%
    filter((IDWE.x==IDWE.y | IDEE.x ==IDEE.y | IDWA.x==IDWA.y | IDEA.x==IDEA.y)) %>%
    filter( !((((!is.na(IDWE.x)) & (!is.na(IDWE.y)) & IDWE.x != IDWE.y)) |
                (((!is.na(IDEE.x)) & (!is.na(IDEE.y)) & IDEE.x != IDEE.y)) |
                (((!is.na(IDWA.x)) & (!is.na(IDWA.y)) & IDWA.x != IDWA.y)) |
                (((!is.na(IDEA.x)) & (!is.na(IDEA.y)) & IDEA.x != IDEA.y)))) %>%
    dplyr::select(n.y) %>%
    unlist() %>%
    as.vector()
  
  
  
  
  
  
  final <- final %>% filter(!n %in% to_remove) 
  
  
  
  
  i=i+1
  
}





final <- final %>% arrange(desc(nb), chr, posl, coloc) %>%
  mutate(clique=seq(1,n(),1)) %>%
  dplyr::select(-n)



duplic <- final %>% dplyr::select(IDEA,IDWA,IDEE,IDWE) %>% pivot_longer(cols = c("IDEA","IDWA","IDEE","IDWE") ) %>%
  dplyr::select(value)%>% na.omit() %>% unlist() %>% as.vector() 

cat("\n\n nb of HR belonging to more than one clique \n\n ")
length(which(duplicated(duplic))) 

# k=length(which(duplicated(duplic)))
# 
# finalf %>% filter(IDWA %in% duplic[which(duplicated(duplic))[k]] | 
#                     IDWE %in% duplic[which(duplicated(duplic))[k]] |
#                     IDEE %in% duplic[which(duplicated(duplic))[k]] |
#                     IDEA %in% duplic[which(duplicated(duplic))[k]] )




single_to_remove <-  final %>% dplyr::select(IDEA,IDWA,IDEE,IDWE) %>% pivot_longer(cols = c("IDEA","IDWA","IDEE","IDWE") ) %>%
  dplyr::select(value)%>% na.omit() %>% unlist() %>% as.vector() %>% unique()

WE <- h %>% filter(population=="WE") %>%
  mutate(clique=NA, nb=1, WE=1, EE=0,WA=0, EA=0, IDWE=IDintHR, IDEE=NA, IDWA=NA, IDEA=NA, coloc="WE") %>%
  dplyr::select(clique, nb, WE, EE, WA, EA, IDWE, IDEE, IDWA, IDEA, posl, posr, chr, coloc)

EE <- h %>% filter(population=="EE") %>%
  mutate(clique=NA, nb=1, WE=0, EE=1,WA=0, EA=0, IDWE=NA, IDEE=IDintHR, IDWA=NA, IDEA=NA, coloc="EE") %>%
  dplyr::select(clique, nb, WE, EE, WA, EA, IDWE, IDEE, IDWA, IDEA, posl, posr, chr, coloc)

WA <- h %>% filter(population=="WA") %>%
  mutate(clique=NA, nb=1, WE=0, EE=0,WA=1, EA=0, IDWE=NA, IDEE=NA, IDWA=IDintHR, IDEA=NA, coloc="WA") %>%
  dplyr::select(clique, nb, WE, EE, WA, EA, IDWE, IDEE, IDWA, IDEA, posl, posr, chr, coloc)

EA <- h %>% filter(population=="EA") %>%
  mutate(clique=NA, nb=1, WE=0, EE=0,WA=0, EA=1, IDWE=NA, IDEE=NA, IDWA=NA, IDEA=IDintHR, coloc="EA") %>%
  dplyr::select(clique, nb, WE, EE, WA, EA, IDWE, IDEE, IDWA, IDEA, posl, posr, chr, coloc)


single <- rbind(WE, EE, WA, EA) %>%
  filter(!(IDWE %in% !!single_to_remove )) %>%
  filter(!(IDEE %in% !!single_to_remove )) %>%
  filter(!(IDWA %in% !!single_to_remove )) %>%
  filter(!(IDEA %in% !!single_to_remove ))


final <- rbind(final, single) %>%
  dplyr::select(clique, coloc, chr, posl, posr, WE,EE,WA,EA, IDWE,IDEE,IDWA,IDEA) %>%
  left_join(h %>% dplyr::select(chr, region, posl), by=c("chr","posl"))  %>%
  dplyr::select(clique, coloc, chr, region,posl, posr, WE,EE,WA,EA, IDWE,IDEE,IDWA,IDEA) %>%
  unique() %>%
  arrange(clique, region) %>%
  group_by(clique, coloc, chr,posl, posr, WE,EE,WA,EA, IDWE,IDEE,IDWA,IDEA) %>%
  slice(1) %>%
  unique() %>%
  ungroup() %>%
  mutate(clique=seq(1,n(),1))




h <- final %>% dplyr::select(IDWE, IDEE, IDWA, IDEA, coloc, clique) %>%
  pivot_longer(cols=c("IDWE","IDEE","IDWA","IDEA"), values_to="IDintHR") %>%
  dplyr::select(coloc, IDintHR, clique) %>% 
  na.omit() %>%
  full_join(h, by="IDintHR", suffix=c("",".x")) %>%
  dplyr::select(population, chr, region, posl, posr, IDintHR, coloc, clique, everything())


################ changer les coloc=NA en coloc=WE

spe <- final %>% group_by(coloc) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(tcliques=sum(n))



vecteur <- c("WEEEWAEA",
             "WEEEWA","WEEEEA","WEWAEA", "EEWAEA",
             "WEEE","WEWA","WEEA","EEWA","EEEA", "WAEA",
             "WE","EE","WA","EA")

if (!identical(sort(as.character(spe$coloc)), sort(vecteur))){
  
  pb <- vecteur[which(!vecteur %in% spe$coloc)]
  
  
  for (i in 1:length(pb)){
    
    
    tab <- data.frame(coloc=pb[i],
                      n=0,
                      tcliques=unique(spe$tcliques))
    
    spe <- rbind(spe, tab)
    
    
    
    
    
    
  }
  
  
  
  
  
}



if (!is.null(unique(fread(titre_HR)$iteration))){
  spe <- spe %>%  mutate(iteration=unique(fread(titre_HR)$iteration))
  
  
}




if  (type == "1000"){
  cat("\n\n dataframe for upset graph (HR_100) \n\n")
  print(head(final))
  write.table(final, titre_sortie_cliques_HR, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
  
} 


if (type=="reference") {
  
  
  cat("\n\n information about HR \n\n")
  print(head(h))
  write.table(h, titre_sortie_cliques_HR, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
  
  
}



cat("\n\n coloc of HR, HR_1000 or permutation \n\n")
spe
write.table(spe, titre_sortie_coloc, col.names = T, row.names = F, dec=".", sep="\t", quote=F)



sessionInfo()




