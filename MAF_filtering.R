


Sys.time()
cat("\n\nMAF_filtering.R\n\n")
rm(list = ls())
set.seed(1)
graphics.off()


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))






variables <- commandArgs(trailingOnly=TRUE)
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_landraces <- variables[1]
titre_WE <- variables[2]
titre_WA <- variables[3]
titre_EE <- variables[4]
titre_EA <- variables[5]
titre_pos <- variables[6]
dossier_sortie <- variables[7]
titre_fid <- variables[8]
titre_SNP_to_keep <- variables[9]
titre_landraces_to_keep <- variables[10]


 # titre_landraces <- "/work/adanguy/these/pipeline/020820/sources/landraces.txt"      
 # titre_WE <- "/work/adanguy/these/pipeline/020820/PLINK/WE.map"               
 # titre_WA <-  "/work/adanguy/these/pipeline/020820/PLINK/EE.map"               
 # titre_EE <- "/work/adanguy/these/pipeline/020820/PLINK/WA.map"               
 # titre_EA <- "/work/adanguy/these/pipeline/020820/PLINK/EA.map"               
 # titre_pos <- "/work/adanguy/these/pipeline/020820/sources/SNP_positions.txt"  
 # dossier_sortie <-  "/work/adanguy/these/pipeline/020820/PLINK/"                     
 # titre_fid <- "/work/adanguy/these/pipeline/020820/PLINK/update_FID.txt"       
 # titre_SNP_to_keep <- "/work/adanguy/these/pipeline/020820/PLINK/SNP_PLINK_3.txt"      
 # titre_landraces_to_keep <- "/work/adanguy/these/pipeline/020820/PLINK/landraces_PLINK_3.txt"
 



# titre_landraces <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/landraces.txt"         
# titre_map <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/Map.map"
# titre_pos <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/SNP_positions.txt"
# dossier_sortie <- "/home/adanguydesd/Documents/These_Alice/pipeline/poubelle/"


cat("\n\n Input 1 : Passeport data of landraces \n\n")
l <- fread(titre_landraces)
head(l)
l <- l %>%
  arrange(LINE) %>%
  filter(kept==T) %>% # remove landraces with > 10% NA and >5% heterozygotie and < 0.5 max admixture and higly related
  mutate(oldFID=1) %>%
  mutate(newFID=population) %>%
  mutate(oldID=LINE) %>%
  mutate(newID=LINE) %>%
  dplyr::select(oldFID, oldID, newFID, newID)


cat("\n\n Input 2 : Physical positions of SNP \n\n")
pos <-fread(titre_pos, header=T) # to have chr and region per SNP
head(pos)
pos <- pos %>%
  mutate(chrregion=paste0(chr, region)) %>%
  dplyr::select(SNP, chrregion)

cat("\n\n Input 3 : SNP in genotyping data \n\n")
mapWE <- fread(titre_WE, header=F) 
head(mapWE)

mapEE <- fread(titre_EE, header=F) 

mapWA <- fread(titre_WA, header=F) 

mapEA <- fread(titre_EA, header=F) 

map <- rbind(mapWE,mapEE,mapWA,mapEA) %>%
  dplyr::select(V1,V2,V4) %>%
  unique() %>%
  rename(chr=V1, SNP=V2, posSNP=V4) %>%
  inner_join(pos, by="SNP") %>%
  arrange(chr, posSNP) %>%
  dplyr::select(SNP, chrregion)

b="3BR1"
for (b in unique(map$chrregion)){
  
  titre_sortie <- paste0(dossier_sortie,"chrregion_", b,".txt")
  liste_SNP_region <- map %>% 
    filter(chrregion==!!b) %>%
    dplyr::select(SNP) %>%
    unlist() %>% 
    as.vector()
  
  if(b==unique(map$chrregion)[1]){
  cat("\n\n Output 1 : SNP per chrregion for PLINK \n\n")
  print(head(liste_SNP_region))
  write.table(liste_SNP_region, titre_sortie, col.names = F, row.names = F, dec=".", sep="\t", quote=F)
  } else {
  write.table(liste_SNP_region, titre_sortie, col.names = F, row.names = F, dec=".", sep="\t", quote=F)
  
  }  
}
  




cat("\n\n Output 1 : update file for PLINK, specifying FID \n\n")
head(l)
write.table(l, titre_fid, col.names = F, row.names = F, dec=".", sep="\t", quote=F)



cat("\n\n Output 2 : final set of SNP \n\n")
map <- map %>% dplyr::select(SNP) %>% unlist() %>% as.vector()
head(map)
write.table(map, titre_SNP_to_keep, col.names = F, row.names = F, dec=".", sep="\t", quote=F)



cat("\n\n Output 3 : final set of landraces \n\n")
l <- l %>% dplyr::select(oldFID, oldID) 
head(l)
write.table(l, titre_landraces_to_keep, col.names = F, row.names = F, dec=".", sep="\t", quote=F)


sessionInfo()