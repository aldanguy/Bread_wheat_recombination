

Sys.time()
cat("\n\nFST_SNP.R\n\n")
rm(list = ls())
set.seed(1)
graphics.off()


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gtools))






variables <- commandArgs(trailingOnly=TRUE)
print(variables)

chrregion <- variables[1]
titre_reynolds <- variables[2]
titre_SNP <- variables[3]
titre_output <- variables[4]


if (chrregion == "1AC"){
  cat("\n\nVariables : \n")
  print(variables)
  cat("\n\n")
}



# chrregion <- "1AC"
# titre_reynolds <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/hapflk_reynolds.txt"
# dossier_SNP <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/"




snp <- fread(titre_SNP, header=F)

if (chrregion == "1AC"){
cat("\n\n Intput 1 : SNP of the chrregion \n\n")
head(snp)
}

nb_SNP <- snp %>% unlist() %>% as.vector() %>% length()

r <- fread(titre_reynolds, header=F, data.table = F)

if (chrregion == "1AC"){
  cat("\n\n Intput 2 : Reynolds distance of the chrregion \n\n")
  head(r)
}

pop <- r$V1
r <- r[,-1]
colnames(r) <- pop
rownames(r) <- pop


noms <- combinations(n = length(pop), r = 2, v = pop, repeats.allowed = TRUE)  
r2 <- data.frame(chrregion, noms, reynold=r[noms]) %>%
  rename(P1=X1, P2=X2) %>%
  mutate(P1.2 = case_when(P1=="WE"|P2=="WE" ~ "WE",
                          (P1!= "WE" & P2 !="WE") & (P1=="EE" | P2 =="EE") ~ "EE",
                          (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE") & (P1=="WA"|P2=="WA") ~ "WA",
                          (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA") & (P2=="EA"|P2=="EA")~ "EA",
                          (P1 != "WE" & P2 !="WE" & P1 !="EE" & P2 != "EE" & P1 !="WA" & P2 != "WA" & P1 !="EA"& P2 !="EA" ) & (P2=="CsRe"|P2=="CsRe")~ "CsRe")) %>%
  mutate(P2.2 = ifelse(P1==P1.2, as.character(P2), as.character(P1))) %>%
  dplyr::select(-one_of("P1","P2")) %>%
  rename(P1=P1.2, P2=P2.2) %>%
  filter(P1 != P2) %>%
  mutate(FST=reynold) %>%
  mutate(methode="SNP") %>%
  mutate(nb_SNP=nb_SNP) %>%
  dplyr::select(chrregion, P1, P2, FST, nb_SNP, methode)







if (chrregion == "1AC"){
  
  cat("\n\n Output 1 : FST of the chrregion computed on SNP \n\n")
  print(head(r2))
  write.table(r2, titre_output, col.names = T, row.names = F, dec=".", sep="\t", quote=F, append = F)
  
  sessionInfo()
  
}else {
  
  write.table(r2, titre_output, col.names = F, row.names = F, dec=".", sep="\t", quote=F, append = T)
  
  
}

sessionInfo()
