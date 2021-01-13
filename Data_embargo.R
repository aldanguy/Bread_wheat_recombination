


Sys.time()
cat("\n\nGraphs_output.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))


titre_280 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/Rimbert_2018.txt"
#titre_map_final <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/landraces_final.map"
titre_map_initial <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/donnees_publiees/BW_261K_wp3.map"
titre_SNP_positions <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/SNP_positions.txt"
# titre_WE <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/WE.map"
# titre_EE <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/EE.map"
# titre_WA <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/WA.map"
# titre_EA <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/EA.map"
# titre_PHASE <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/PHASE_summary_outputs.txt"
# 
titre_snp_pos_280k
titre_mqs_280k_initial <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/SNP_PLINK4_mqs_280k_initial"

# proportion of SNP in final map from TABW280k

liste_280k <- fread(titre_280) %>% dplyr::select(PROBESET_ID) %>% unlist() %>% as.vector()
snp_pos <- fread(titre_SNP_positions)

# map_final <- fread(titre_map_final)
map_initial <- fread(titre_map_initial)
# WE<- fread(titre_WE)
# EE<- fread(titre_EE)
# WA<- fread(titre_WA)
# EA <- fread(titre_EA)
# PHASE <- fread(titre_PHASE) %>% filter(w_center==T)
# 
# length(which(map_final$V2 %in% liste_280k))/nrow(map_final)
# length(which(WE$V2 %in% liste_280k))/nrow(WE)
# length(which(EE$V2 %in% liste_280k))/nrow(EE)
# length(which(WA$V2 %in% liste_280k))/nrow(WA)
# length(which(EA$V2 %in% liste_280k))/nrow(EA)
# 
# 
# 
# PHASE %>% filter(population=="WE") %>% filter(SNPlpop %in% liste_280k) %>% nrow()/PHASE %>% filter(population=="WE") %>% nrow()
# PHASE %>% filter(population=="EE") %>% filter(SNPlpop %in% liste_280k & SNPrpop %in% liste_280k) %>% nrow()/PHASE %>% filter(population=="EE") %>% nrow()
# PHASE %>% filter(population=="WA") %>% filter(SNPlpop %in% liste_280k & SNPrpop %in% liste_280k) %>% nrow()/PHASE %>% filter(population=="WA") %>% nrow()
# PHASE %>% filter(population=="EA") %>% filter(SNPlpop %in% liste_280k & SNPrpop %in% liste_280k) %>% nrow()/PHASE %>% filter(population=="EA") %>% nrow()
# 


liste_mqs_280k_initial <- map_initial$V2[which(map_initial$V2 %in% liste_280k)]
snp_pos_280k <- snp_pos %>% filter(snp_pos$SNP %in% liste_mqs_280k_initial)

write.table(liste_mqs_280k_initial, titre_mqs_280k_initial, col.names = F, row.names=F, dec=".", sep="\t", quote=F)
write.table(snp_pos_280k, titre_snp_pos_280k, col.names = T, row.names=F, dec=".", sep="\t", quote=F)

sessionInfo()