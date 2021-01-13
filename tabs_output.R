
#### Formatting of tables

Sys.time()
cat("\n\nTables_output.R\n\n")
rm(list = ls())
graphics.off()
set.seed(1)


suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))



variables <- commandArgs(trailingOnly=TRUE)

cat("\n Graphics and stats\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")

repertoire <- variables[1]
titre_resume <-variables[2]
titre_PHASE_summary_outputs2 <- variables[3]
titre_csre_genetic_map <- variables[4]
titre_landraces_genetic_map <- variables[5]
titre_genetic_maps <- variables[6]
titre_hr <- variables[7]
titre_hr2 <- variables[8]
titre_cliques <- variables[9]
titre_cliques2 <- variables[10]
titre_rho <- variables[11]
titre_lambda <- variables[12]
titre_correlations_mixed_models <- variables[13]
titre_landraces_4Mb <- variables[14]
titre_csre_4Mb <- variables[15]
titre_cor_4Mb <- variables[16]
titre_maps_4Mb <- variables[17]
titre_permutations_HR <- variables[18]
titre_overlap_HR <- variables[19]
titre_overlap_HR_published <-variables[20]

# titre_csre_genotyping <- variables[18]
# titre_csre_genotyping_matrix <- variables[19]
titre_FST_SNP <- variables[21]
titre_FST_haplotypic_blocs <- variables[22]
titre_FST <- variables[23]
titre_SNP_positions <- variables[24]
titre_chr_regions <- variables[25]
titre_correspondance_chr <- variables[26]
titre_stab1 <- variables[27]
titre_landraces <- variables[28]
titre_stab2 <- variables[29]

# titre_resume <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/PHASE_summary_outputs.txt"
# titre_csre_genetic_map <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/csre_genetic_map.txt"
# titre_landraces_genetic_map <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/landraces_genetic_map.txt"
# titre_FST_SNP <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/chrregion_FST_SNP.txt"
# titre_FST_haplotypic_blocs <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/chrregion_FST_haplotypic_blocks.txt"
# titre_SNP_positions <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/SNP_positions.txt"
# titre_chr_regions <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/Decoupage_chr_ble.tab"
# titre_correspondance_chr <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/Codes_chr.txt"
# titre_hr <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/HR.txt"
# titre_cliques <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/cliques_HR.txt"
# titre_lambda <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/chrregion_lambda_random.txt"
# titre_lambda2 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/common_SNP/chrregion_lambda_random.txt"
# titre_rho <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/chrregion_rho_random.txt"
# titre_rho2 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/common_SNP/chrregion_rho_random.txt"
# titre_correlations_mixed_models2 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/common_SNP/correlations_mixed_models_published.txt"
# titre_landraces_4Mb <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/map_4mb_landraces.txt"
# titre_csre_4Mb <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/map_4mb_csre.txt"
# titre_cor_4Mb<- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/correlations_4mb_published.txt"
# titre_landraces <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/landraces.txt"
# titre_csre_genotyping <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/csre_genotyping.txt"
# 
# titre_overlap_HR <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/coloc_HR.txt"
# 
# titre_permutations_HR <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/permutations.txt"
# titre_overlap_HR_published <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/tabs/overlap_HRIs_published.txt"

head(fread(titre_resume))
head(fread(titre_csre_genetic_map))
head(fread(titre_landraces_genetic_map))
head(fread(titre_hr))
head(fread(titre_cliques))
head(fread(titre_lambda))
head(fread(titre_rho))
head(fread(titre_landraces_4Mb))
head(fread(titre_csre_4Mb))
head(fread(titre_permutations_HR))
head(fread(titre_overlap_HR))
head(fread(titre_FST_SNP))



map_4Mb <- fread(titre_landraces_4Mb) %>%
  dplyr::select(-sd_lambda_rho) %>%
  rename(rec_rate=lambda_rho) %>%
  rbind(.,fread(titre_csre_4Mb) %>% dplyr::select(-sdCbay) %>% rename(rec_rate=cbay)) %>%
  mutate(population=factor(population, levels=c("CsRe","WE","EE","WA","EA"))) %>%
  arrange(population, chr, posl)
  




# PHASE summary outputs
res <- fread(titre_resume) %>% 
  filter(w_center==T) %>%
  dplyr::select(population, chr, region, ID, posSNPlpop, posSNPrpop, SNPlpop, SNPrpop, lambda_med, rho_med, lambda_rho_med) %>%
  rename(PHASE_window_ID=ID,
         SNPl=SNPlpop,
         SNPr=SNPrpop,
         posl=posSNPlpop,
         posr=posSNPrpop) %>%
  arrange(population, chr, posl)



### genetic maps



genetic_maps <- fread(titre_csre_genetic_map) %>% 
  dplyr::select(population, chr, region, posSNPlpop, posSNPrpop, SNPlpop, SNPrpop, cbay) %>%
  rename( SNPl=SNPlpop,
          SNPr=SNPrpop,
          posl=posSNPlpop,
          posr=posSNPrpop,
          meiotic_rec_rate=cbay) %>%
  mutate(historical_rec_rate=NA, coef_proportionnality=NA) %>%
  dplyr::select(population, chr, region, posl, posr, SNPl, SNPr, historical_rec_rate, meiotic_rec_rate, coef_proportionnality) %>%
  rbind(.,fread(titre_landraces_genetic_map)) %>%
  dplyr::select(population, chr, region, posl, posr, SNPl, SNPr, meiotic_rec_rate, historical_rec_rate, coef_proportionnality) %>%
  arrange(population, chr, posl)




# hr

hr <- fread(titre_hr) %>%
  filter(kept==T) %>%
  dplyr::select(population, chr, region, posl, posr, IDintHR, lambda_med) %>%
  rename(HRI_ID=IDintHR) %>%
  arrange(population, posl)



# cliques
cliques <- fread(titre_cliques) %>%
  dplyr::select(chr, region, posl, posr, coloc, clique) %>%
  group_by(chr, region, coloc, clique) %>%
  summarise(posl=max(posl), posr=min(posr)) %>%
  rename(populations=coloc) %>%
  dplyr::select(chr, region, posl, posr, populations) %>%
  arrange(chr, posl)   %>%
  ungroup() %>%
  mutate(ID_clique=paste0("clique_",1:nrow(.)))



# mixed model correlations

lambda <- fread(titre_lambda) %>%
  dplyr::select(chrregion, P1, P2, estimateur, BIC) %>%
  mutate(chr=substr(chrregion, 1,2)) %>%
  mutate(region=substr(chrregion, 3,5)) %>%
  rename(cor_lambda=estimateur,
         BIC_lambda=BIC) %>%
  dplyr::select(chr, region, P1, P2, cor_lambda, BIC_lambda)


correlations_mixed_models <- fread(titre_rho) %>%
  dplyr::select(chrregion, P1, P2, estimateur, BIC) %>%
  mutate(chr=substr(chrregion, 1,2)) %>%
  mutate(region=substr(chrregion, 3,5)) %>%
  rename(cor_rho=estimateur,
         BIC_rho=BIC) %>%
  dplyr::select(chr, region, P1, P2, cor_rho, BIC_rho) %>%
  full_join(lambda, by=c("chr"="chr","region"="region", "P1"="P1", "P2"="P2")) %>%
  dplyr::select(chr, region, P1, P2, cor_lambda, BIC_lambda, cor_rho, BIC_rho)  %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
  mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
  arrange(chr, region, P1, P2) %>%
  filter(P1 !=P2)






# 4Mb correlation 


cor_4Mb<- fread(titre_landraces_4Mb) %>%
  inner_join(fread(titre_csre_4Mb), by=c("chr"="chr","posl"="posl", "posr"="posr", "region"="region")) %>%
  rename(population=population.x) %>%
  group_by(population, chr, region) %>%
  summarise(correlation=cor(lambda_rho, cbay), nb_windows_4Mb=n())  %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  arrange(population, chr, region)


# overlap of HR

overlap_HR <- fread(titre_permutations_HR) %>%
  rename(populations=coloc, nbHRIs=n, nb_total_cliques=tcliques) %>%
  mutate(iteration=iteration+1) %>%
  rbind(fread(titre_overlap_HR) %>%
          rename(populations=coloc, nbHRIs=n, nb_total_cliques=tcliques) %>%
          mutate(iteration="ref"))  %>%
  arrange(iteration, populations) 



# Fst
FST_SNP <- fread(titre_FST_SNP) %>% dplyr::select(chrregion, P1, P2, FST, methode)
FST <- fread(titre_FST_haplotypic_blocs) %>% dplyr::select(chrregion, P1, P2, FST, methode) %>%
  rbind(., FST_SNP) %>%
  mutate(chr=substr(chrregion, 1,2)) %>%
  mutate(region=substr(chrregion, 3,5)) %>%
  rename(method=methode) %>%
  mutate(method=ifelse(method=="haplotype", "haplotypic_blocks", method)) %>%
  dplyr::select(chr, region, method, P1, P2, FST) %>%
  mutate(region=factor(region, levels=c("R1","R2a","C","R2b","R3"))) %>%
  mutate(P1=factor(P1, levels=c("WE","EE","WA","EA"))) %>%
  mutate(P2=factor(P2, levels=c("WE","EE","WA","EA"))) %>%
  arrange(chr, region, method, P1, P2)





if (repertoire=="all"){
  
  head(fread(titre_landraces))
  head(fread(titre_SNP_positions))
  head(read.table(titre_correspondance_chr, header=F, dec=".", sep="\t"))
  head(read.table(titre_chr_regions, header=T, dec=".", sep="\t", skip=1))

  
  
  # landraces data
  
  
  

  
  landraces <- fread(titre_landraces) %>%
    filter(kept==T) %>%
    dplyr::select(population, LINE, Nom_TaBW420K, country, STRgp4, STRarea8, STRmemb4_3, STRmemb4_4,STRmemb4_1,STRmemb4_2) %>%
    rename(ID_TaBW420K=Nom_TaBW420K,
           ID_LINE=LINE,
           STRgp8=STRarea8,
           CAA=STRmemb4_1,
           SEA=STRmemb4_2,
           NWE=STRmemb4_3,
           SEE=STRmemb4_4) %>%
    mutate(STRgp4=case_when(STRgp4==2 ~ "CAA",
                            STRgp4==5 ~ "SEA",
                            STRgp4==8 ~ "NWE",
                            STRgp4==12 ~ "SEE"))  %>%
    mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
    arrange(population, ID_LINE)
  
  

  
  
  # # csre genotyping
  # geno <- fread(titre_csre_genotyping) %>%
  #   dplyr::select(-population, -region) %>%
  #   rename(pos=posSNP) %>%
  #   arrange(chr, pos)
  
  
  ### genomic regions
  

  
  
  
  cor_chr2 <- read.table(titre_correspondance_chr, header=F, dec=".", sep="\t") %>%
    rename( chr = V1, chr_code_nombre = V2) %>%
    mutate(chr=as.character(chr))
  head(cor_chr2)
  
  
  regions2  <- read.table(titre_chr_regions, header=T, dec=".", sep="\t", skip=1) %>%
    dplyr::select(Chromosome, R1.R2a, R2a.C, C.R2b, R2b.R3) %>%
    mutate(chr = str_remove(Chromosome,"chr")) %>%
    dplyr::select(-one_of("Chromosome")) %>% 
    full_join(cor_chr2, by="chr") %>%
    rename(R1=R1.R2a) %>%
    rename(R2a=R2a.C) %>%
    rename(C=C.R2b) %>%
    rename(R2b=R2b.R3) %>%
    mutate(R3=1e9) %>% # artificial end of R3
    mutate(R1=R1*1e6) %>%
    mutate(R2a=R2a*1e6) %>%
    mutate(C=C*1e6) %>%
    mutate(R2b=R2b*1e6) %>%
    pivot_longer(-c("chr_code_nombre","chr"), names_to ="region", values_to = "pos_max") %>%
    group_by(chr_code_nombre) %>%
    mutate(pos_min=lag(pos_max) +1) %>%
    mutate(pos_min=ifelse(!is.na(pos_min),pos_min,1)) %>%
    mutate(chrregion=paste0(chr,region)) %>%
    mutate(pos=pos_min+((pos_max-pos_min)/2)) %>%
    ungroup() %>%
    dplyr::select(chr, region, pos_min, pos_max) %>%
    rename(posl=pos_min, posr=pos_max)
  
  head(regions2)
  
  
  genomic_regions <- fread(titre_SNP_positions) %>%
    filter(region=="R3") %>%
    group_by(chr, region) %>%
    summarise(posmax=max(posSNP)) %>% 
    full_join(regions2, by=c("chr"="chr","region"="region")) %>%
    mutate(posr=ifelse(region=="R3",posmax,posr)) %>%
    mutate(size_Mb=round((posr-posl)/1e6)) %>%
    dplyr::select(chr, region, posl, posr, size_Mb) %>%
    arrange(chr, posl)
  
  
  
  print(head(genomic_regions))
  write.table(genomic_regions, titre_stab1, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
  # geno[1:10,1:10]
  # write.table(geno, titre_csre_genotyping_matrix, col.names = T, row.names = F, dec=".", sep="\t", quote=F)

  print(head(landraces))
  write.table(landraces, titre_stab2, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
  
  
  
  
  
  
  
}



head(cor_4Mb)
write.table(cor_4Mb, titre_cor_4Mb, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
head(correlations_mixed_models)
write.table(correlations_mixed_models, titre_correlations_mixed_models, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
head(cliques)
write.table(cliques, titre_cliques2, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
head(hr)
write.table(hr, titre_hr2, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
head(genetic_maps)
write.table(genetic_maps, titre_genetic_maps, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
# column coef_proportionnality (no unit) = average ratio bewteen historical_rec_rate (converted from 4Ne*M/pb to 4*Ne cM/Mb , i.e. multiplied by a 1e-8 factor) divided by average Csre Bayesian meiotic rate in the genomic region (1AR1...7DR3)
# column meiotic_rec_rate (cM/Mb) : for CsRe = Bayesian meiotic rec rate ; for WE, EE, WA and EA =  historical_rec_rate * 1e8 / coef_proportionnality
head(res)
write.table(res, titre_PHASE_summary_outputs2, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
# Outputs from PHASE, median per interval
head(map_4Mb)
write.table(map_4Mb, titre_maps_4Mb, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
# Outputs from PHASE, median per interval
head(overlap_HR)
write.table(overlap_HR, titre_overlap_HR_published, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
head(FST)
write.table(FST, titre_FST, col.names = T, row.names = F, dec=".", sep="\t", quote=F)

sessionInfo()