Sys.time()
cat("\n\nPHASE_2.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))



variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_monitor <- variables[1]
# Contient toutes les sorties _monitor de PHASE

titre_fenetres <- variables[2]
# Contient des infos sur les fenetres rentrées dans PHASE
# Produit par Fenetres_groupe3_140619.R

titre_scaffolds <- variables[3]
titre_phase <- variables[4]
titre_pos <- variables[5]
titre_sortie <- variables[6]
titre_mqs <- variables[7]
titre_resume <- variables[8]
titre_png <- variables[9]




# titre_monitor <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/WE_3B_28_monitor"
# titre_fenetres <- "/home/adanguydesd/Documents/These_Alice/pipeline/compa_orge/180619/sorties_PHASE/bilan/synthese_fenetres_groupe3.txt"
# titre_scaffolds <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/scaffolds.txt"
# titre_phase <-"/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/WE_3B_28_sorties_brutes_PHASE.txt"
# titre_pos <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/SNP_positions.txt"
# titre_sortie <- "/home/adanguydesd/Documents/These_Alice/pipeline/fichiers_test/sorties_PHASE_regroupees/cluster_1_1_3A_ID100ID_sorties_PHASE_traitees.txt"

# titre_monitor <-"/work/adanguy/these/pipeline/020820/PHASE/output/raw/WA_1A_53_monitor"                      
# titre_fenetres <- "/work/adanguy/these/pipeline/020820/PHASE/windows/description/w_WE_1A.txt"                 
# titre_scaffolds <- "/work/adanguy/these/pipeline/020820/sources/scaffolds.txt"                                 
# titre_phase <- "/work/adanguy/these/pipeline/020820/PHASE/output/format/WE_1A_9_PHASE_raw_outputs.txt"     
# titre_pos <-  "/work/adanguy/these/pipeline/020820/sources/SNP_positions.txt"                             
# titre_sortie <- "/work/adanguy/these/pipeline/020820/PHASE/output/format/WE_1A_9_PHASE_raw_outputs.txt"     
# titre_mqs <- "/work/adanguy/these/pipeline/020820/intervals/mqs/WE_1A_9_mqs.txt"                         
# titre_resume <-"/work/adanguy/these/pipeline/020820/PHASE/output/summary/WE_1A_9_PHASE_summary_outputs.txt"
# titre_png <-"/work/adanguy/these/pipeline/020820/graphs/PHASE/WE_1A_9.png"  

cat("\n\n Input 1 : Historical recombination estimates for the window \n\n")
phase <- fread(titre_phase, header=T, dec=".", sep="\t", data.table = F) 
head(phase)
phase <- phase %>%
  dplyr::select(population, chr, ID, posSNPlpop, posSNPrpop, rho, lambda, nb_underflows, nb_fails) %>%
  mutate(intervalle=paste0(chr, "_", posSNPlpop, "_", posSNPrpop)) %>%
  group_by(intervalle) %>%
  mutate(num_row_PHASE_output = seq(1,min(table(intervalle)),1)) 

posterior <- max(phase$num_row_PHASE_output)


##################################### A changer
# Pour savoir les l'intervalle est inclue dans la compa
# fen <- fread(titre_fenetres, header=T, dec=".", sep="\t", data.table = F) %>%
#   dplyr::select(chr, centerlfen, centerrfen, ID) %>%
#   filter(ID == unique(phase$ID))
cat("\n\n Input 2 : Windows limits \n\n")
fen <- fread(titre_fenetres, header=T, dec=".", sep="\t", data.table = F) 
head(fen)
fen <- fen %>%
   dplyr::select(chr, centerlfen, centerrfen, ID) %>% ##################################################################
   filter(ID == unique(phase$ID))







# Pour savoir si l'intervalle est entre deux scaffolds
cat("\n\n Input 3 : scaffolds limits \n\n")
scaffolds <- fread(titre_scaffolds, header=T, dec=".", sep="\t", data.table = F) 
head(scaffolds)
scaffolds <- scaffolds %>%
  dplyr::select(chr, poslscaf, posrscaf) %>%
  filter(chr == unique(phase$chr))

cat("\n\n Input 4 : physical positions of SNP \n\n")
pos <- fread(titre_pos, header=T, dec=".", sep="\t", data.table = F) 
head(pos)
pos <- pos %>%
  dplyr::select(SNP,chr,region,posSNP) %>%
  filter(chr == unique(phase$chr))

cat("\n\n Input 5 : likelihood variation during the main iterations \n\n")
monitor <- fread(titre_monitor, header=F, dec=".", sep=" ", data.table = F)
head(monitor)

#Certains fichiers peuevnt être vides. Dans ce cas, ne pas les analyser
if (min(monitor$V2) > -Inf){
  
  ID=strsplit(titre_monitor, split="/") %>%
    unlist() %>%
    last() %>%
    gsub("_monitor","",.)
  
  
  
  monitor <-  monitor %>% mutate(index=seq(1,nrow(.), 1)) 
  
  modele <- lm(V2~index, data=monitor)
  
  PAC_slope <- as.vector(modele$coefficients[2])
  pente=round(PAC_slope, digits=2)
  
  cat("\n\n Graph 1 : likelihood variation during the main iterations \n\n")
  png(titre_png)
  
  plot(monitor$index, monitor$V2, main= ID, type="b", xlab="iteration", ylab="PAC-B likelihood", pch=19)
  abline(modele$coefficients, col="red", lwd=2)
  legend("bottom", legend=pente, col="red", lty=1)
  
  dev.off()
  
} else {
  
  PAC_slope <- NA
  
  
}

phase[,"PAC_slope"] <- PAC_slope


different_scaffold <- phase %>%
  left_join(scaffolds, by="chr") %>%
  filter(!( posSNPlpop >= posrscaf | posSNPrpop <= poslscaf))  %>%
  group_by(chr, posSNPlpop, posSNPrpop) %>%
  summarise(nscaf=n()/posterior) %>%  # catch intervals whose SNP are located on two different scaffolds
  unique() %>%
  mutate(different_scaffold=ifelse(nscaf==1, F, T)) %>%
  dplyr::select(chr, posSNPlpop, posSNPrpop, different_scaffold)



phase <- phase %>%
  full_join(different_scaffold, by=c("chr","posSNPlpop","posSNPrpop")) %>%
  mutate(w_center=ifelse(posSNPlpop>=fen$centerlfen & posSNPrpop <= fen$centerrfen, "TRUE","FALSE")) %>%
  left_join(pos, by=c("posSNPlpop"="posSNP"), suffix = c("", ".x")) %>%
  left_join(pos, by=c("posSNPrpop"="posSNP"), suffix = c("", ".y")) %>%
  mutate(methode="landraces") %>%
  mutate(lambda_rho = rho*lambda) %>%
  rename(SNPlpop=SNP, SNPrpop=SNP.y) %>%
  arrange(posSNPlpop, num_row_PHASE_output) %>%
  mutate(lpop=posSNPrpop-posSNPlpop) %>%
  mutate(lpop=lpop/1e6) %>%
  ungroup() %>%
  dplyr::select(population,
                chr,
                region,
                ID,
                SNPlpop,
                SNPrpop,
                posSNPlpop,
                posSNPrpop,
                lpop,
                different_scaffold,
                lambda,
                rho,
                lambda_rho,
                num_row_PHASE_output,
                w_center,
                PAC_slope,
                nb_fails,
                nb_underflows,
                methode) 

phase <- phase %>% mutate(lambda=ifelse((length(unique(phase$lambda))==length(unique(phase$SNPlpop)) | is.na(PAC_slope)), NA, lambda)) %>%
  mutate(lambda_rho=ifelse((length(unique(phase$lambda))==1 | length(unique(phase$rho))==1 | is.na(PAC_slope)), NA, lambda_rho)) %>%
  mutate(rho=ifelse((length(unique(phase$rho))==1 | is.na(PAC_slope)), NA, rho))  %>% as.data.frame()




cat("\n\n Output 1 : PHASE outputs formated and with additionnal info  \n\n")
head(phase)
write.table(phase,titre_sortie, col.names = T, dec=".",sep="\t", row.names = F, quote = F)






rho_med <- median(phase$rho)
rho_mean <- mean(phase$rho)
rho_sd <- sd(phase$rho)


resume <- phase %>% 
  group_by(population,
           chr,
           region,
           ID,
           SNPlpop,
           SNPrpop,
           posSNPlpop,
           posSNPrpop,
           lpop,
           different_scaffold,
           w_center,
           PAC_slope,
           nb_fails,
           nb_underflows,
           methode) %>%
  summarise(lambda_rho_med = median(lambda_rho),
         lambda_rho_mean = mean(lambda_rho),
         lambda_mean = mean(lambda),
         lambda_med=median(lambda),
         lambda_sd=sd(lambda),
         lambda_rho_sd=sd(lambda_rho)) %>%
  mutate(rho_med = rho_med,
         rho_mean=rho_mean,
         rho_sd=rho_sd) %>%
  arrange(posSNPlpop) %>%
  dplyr::select(population,
                chr,
                region,
                ID,
                SNPlpop,
                SNPrpop,
                posSNPlpop,
                posSNPrpop,
                lpop,
                different_scaffold,
                lambda_med,
                lambda_mean,
                lambda_sd,
                rho_med,
                rho_mean,
                rho_sd,
                lambda_rho_med,
                lambda_rho_mean,
                lambda_rho_sd,
                w_center,
                PAC_slope,
                nb_fails,
                nb_underflows,
                methode) %>%
  unique() %>%
  as.data.frame()


cat("\n\n Output 2 : PHASE outputs summarized and with additionnal info  \n\n")
head(resume) 
write.table(resume, titre_resume, col.names = T, dec=".",sep="\t", row.names = F, quote = F)


# Can be empty at three conditions : PAC_slope is NA OR all lambda are equal to 1 or all rho are the same
mqs <- phase %>%
  filter(w_center==T) %>%
  filter(complete.cases(.)) %>%
  dplyr::select(SNPlpop, SNPrpop) %>%
  gather() %>%
  dplyr::select(value) %>%
  unlist() %>%
  as.vector() %>%
  unique() 
if (length(mqs) >0){
  
  cat("\n\n Output 3 : marqueurs to define de-novo intervals  \n\n")
  print(head(mqs))
  write.table(mqs, titre_mqs, col.names = F, dec=".",sep="\t", row.names = F, quote = F)
  
  
}







cat("\n\n")
sessionInfo()
cat("\n\n")
