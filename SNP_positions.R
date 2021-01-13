
Sys.time()
cat("\n\nSNP_positions.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))





variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")



titre_pos <-variables[1]
titre_SNP_category <-variables[2]
titre_correspondance_chr <-variables[3]
titre_chr_tab <-variables[4]
titre_autres_SNP <- variables[5]
titre_BW <- variables[6]
titre_sortie <- variables[7]


 
     # titre_pos <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Vraies_positions_marqueurs.txt"
     # titre_SNP_category <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Axiom_TaBW420_OrigineSNPs.csv"
     # titre_correspondance_chr <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Codes_chr.txt"
     # titre_chr_tab <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/Decoupage_chr_ble.tab"
     # titre_autres_SNP <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/CartPolyHigh-Oriented_def_def_def.txt"
     # titre_BW <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/landraces.map"
     # titre_sortie <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/SNP_positions.txt"




titre_SNP_category <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/Axiom_TaBW420_OrigineSNPs.csv"
titre_pos <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/Vraies_positions_marqueurs.txt"
titre_correspondance_chr <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/Codes_chr.txt"
titre_autres_SNP <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/CartPolyHigh-Oriented_def_def_def.txt"
titre_BW <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/donnees_publiees/BW_261K_wp3.map"
titre_chr_tab <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/Decoupage_chr_ble.tab"

 # titre_pos <-"/work/adanguy/these/pipeline/amont/Vraies_positions_marqueurs.txt"       
 # titre_SNP_category <-"/work/adanguy/these/pipeline/amont/Axiom_TaBW420_OrigineSNPs.csv"        
 # titre_correspondance_chr <- "/work/adanguy/these/pipeline/amont/Codes_chr.txt"                        
 # titre_chr_tab <- "/work/adanguy/these/pipeline/amont/Decoupage_chr_ble.tab"                
 # titre_autres_SNP <-"/work/adanguy/these/pipeline/amont/CartPolyHigh-Oriented_def_def_def.txt"
 # titre_BW <- "/work/adanguy/these/pipeline/amont/BW_261K_wp3.map"                      
 # titre_sortie <- "/work/adanguy/these/pipeline/160320/sources/SNP_positions.txt" 
 


chr_tab  <- read.table(titre_chr_tab, header=T, dec=".", sep="\t", skip=1) %>%
  dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
  mutate(chr = str_remove(Chromosome,"chr")) %>%
  dplyr::select(-one_of("Chromosome"))
  


correspondance_chr <- read.table(titre_correspondance_chr, header=F, dec=".", sep="\t") %>%
  rename( chr = V1, chr_code_nombre = V2)


mqs_publiques <- fread(titre_SNP_category, header=T, dec=".", sep="\t", data.table = F) %>% 
  rename(SNP = PROBESET_ID, other_SNP_ID = SNP_ID, category = ORIGIN) %>%
  filter(! category %in% c("BGA_validated", "BW_genes"))


autres_SNP <- fread(titre_autres_SNP, header=T, dec=".", sep="\t", data.table = F) %>%
  filter(!is.na(BW_published)) %>%
  rename(pos_genetique=pos, chr1=chr)%>%
  rename(chr=chr_published, chr_code_nombre=chr2_published, posSNP=pos_published, SNP=BW_published) %>%
  dplyr::select(one_of("SNP", "chr", "chr_code_nombre", "posSNP"))


BW <- fread(titre_BW, header=F, dec=".", sep="\t", data.table = F) %>%
  dplyr::select(V2) %>%
  unlist() %>%
  as.vector()


# Importation des vraies positions des mqs
# 1) supprimer les lignes où il y a un NA
# 2) vérifier que chr et chr_coe_nombre designent le mm chr
# 3) supprimer les mqs positionnes aux memes endroits
pos <- fread(titre_pos, header=T, dec=".", sep="\t", data.table = F) %>%
  na.omit() %>%
  rename(SNP=BW, chr=chr, chr_code_nombre=chr2, posSNP=pos)


pos1 <- rbind(pos, autres_SNP) %>%
  unique() %>%
  inner_join(correspondance_chr, by="chr_code_nombre", suffix = c("", "_check")) %>%
  filter(chr == chr_check) %>%
  mutate(ID_pos = paste0(chr, posSNP)) %>%
  filter(! ID_pos %in% unique(.[["ID_pos"]][duplicated(.[["ID_pos"]])])) %>%
  inner_join(mqs_publiques, by="SNP", suffix = c("", ".")) %>%
  inner_join(chr_tab, by="chr", suffix = c("", ".")) %>%
  mutate(region = ifelse(posSNP <=  R1.R2a*1e6, "R1",
                         ifelse(posSNP >= R1.R2a*1e6 +1 & posSNP <= R2a.C*1e6, "R2a",
                                ifelse(posSNP >= R2a.C*1e6 +1 & posSNP <= C.R2b*1e6, "C", 
                                       ifelse(posSNP >= C.R2b*1e6 +1 & posSNP <= R2b.R3*1e6, "R2b",
                                              ifelse(posSNP >= R2b.R3*1e6 +1, "R3", NA)))))) %>%
  dplyr::select(chr,chr_code_nombre, region, SNP, posSNP,other_SNP_ID,category) %>%
  filter(SNP %in% !!BW | SNP %in% autres_SNP$SNP) %>%
  arrange(chr, posSNP)

  

head(pos1)

# which(!pos1$SNP %in% mqs_publiques$SNP)
# 
# which(!fread(titre_autres_SNP, header=T, dec=".", sep="\t", data.table = F) %>%
#   filter(!is.na(BW_published)) %>%
#   dplyr::select(BW_published) %>%
#   unlist() %>%
#   as.vector() %in% mqs_publiques$SNP)
# 
# 
# length(grep("OTV",BW[which(!BW %in% mqs_publiques$SNP)]))
# length(BW[which(!BW %in% mqs_publiques$SNP)])



# pos2 <- rbind(pos, autres_SNP) %>%
#   unique() %>%
#   inner_join(correspondance_chr, by="chr_code_nombre", suffix = c("", "_check")) %>%
#   filter(chr == chr_check) %>%
#   mutate(ID_pos = paste0(chr, posSNP)) %>%
#   inner_join(mqs_publiques, by="SNP", suffix = c("", ".")) %>%
#   inner_join(chr_tab, by="chr", suffix = c("", ".")) %>%
#   mutate(region = ifelse(posSNP <=  R1.R2a*1e6, "R1",
#                          ifelse(posSNP >= R1.R2a*1e6 +1 & posSNP <= R2a.C*1e6, "R2a",
#                                 ifelse(posSNP >= R2a.C*1e6 +1 & posSNP <= C.R2b*1e6, "C", 
#                                        ifelse(posSNP >= C.R2b*1e6 +1 & posSNP <= R2b.R3*1e6, "R2b",
#                                               ifelse(posSNP >= R2b.R3*1e6 +1, "R3", NA)))))) %>%
#   dplyr::select(chr,chr_code_nombre, region, SNP, posSNP,other_SNP_ID,category) %>%
#   arrange(chr, posSNP)

  







##### Nettoyage matrice : garder tous les mqs de "published" et tous les mqs publiques
# Fichier régions
# Fichier intervalles CsRe à l'image de  dossier_intervalles <- "/home/adanguydesd/Documents/These_Alice/pipeline/compa_orge/160719/fichiers_sources/"
# Fichier "update"
# Fichier PLINK

# Autre script pour calculer distance genetique CsRe
# Autre script pour calculer interpolation
# Autre script pour calculer carte CsRe bayésien
  
  
  
########### V2



# titre_pos <- "/home/adanguydesd/Documents/These_Alice/Fichiers_ref_ble/marqueurs/Fichiers_sources/Axiom_TaBW420_Physical_Position_refseqv2.txt"
#  
#  pos <- fread(titre_pos, header=T, dec=".", sep="\t", data.table = F) %>%
#    rename(other_SNP_ID = SNP_ID, chr_code_lettre = chr, chr_code_nombre = chr_code, pos=position) %>%
#    filter(! chr_code_lettre=="uk") %>%
#    mutate(chr_code_nombre = as.numeric(as.character(chr_code_nombre))) %>%
#    inner_join(correspondance_chr, by="chr_code_nombre", suffix = c("", "_check")) %>%
#    filter(chr_code_lettre == chr_code_lettre_check) %>%
#    mutate(ID_pos = paste0(chr_code_lettre, pos)) %>%
#    filter(! ID_pos %in% unique(.[["ID_pos"]][duplicated(.[["ID_pos"]])])) %>%
#    inner_join(mqs_publiques, by="SNP", suffix = c("", ".")) %>%
#    inner_join(chr_tab, by="chr_code_lettre", suffix = c("", ".")) %>%
#    mutate(region = ifelse(pos <=  R1.R2a*1e6, "R1",
#                           ifelse(pos >= R1.R2a*1e6 +1 & pos <= R2a.C*1e6, "R2a",
#                                  ifelse(pos >= R2a.C*1e6 +1 & pos <= C.R2b*1e6, "C", 
#                                         ifelse(pos >= C.R2b*1e6 +1 & pos <= R2b.R3*1e6, "R2b",
#                                                ifelse(pos >= R2b.R3*1e6 +1, "R3", NA)))))) %>%
#    dplyr::select(one_of("SNP", "chr_code_lettre", "region", "pos","chr_code_nombre", "other_SNP_ID", "ctx_strand","ctx_start","ctx_stop", "ctx_coverage", "ctx_identity", "ctx_snp_pos", "category")) %>%
#    arrange(chr_code_lettre, pos)
#  
#  
#  


write.table(pos1, titre_sortie, sep="\t", col.names = T, row.names = F, dec=".", quote=F)
#write.table(pos2, titre_sortie_tabw, sep="\t", col.names = T, row.names = F, dec=".", quote=F)


sessionInfo()


