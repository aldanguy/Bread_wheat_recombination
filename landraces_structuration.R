
Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)




cat("\n\nlandraces_structuration\n\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_landraces <- variables[1]
titre_dissimilarities <- variables[2]
titre_landraces_to_keep <- variables[3]
dossier_sortie_liste <- variables[4]
dossier_graphes <- variables[5]






## Rq : interaction negative entre data.table et dendextend

suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qgraph))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(qvalue, lib.loc = "/work/adanguy/R_library", quietly = T))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gplots))
#suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gridGraphics))
suppressPackageStartupMessages(library(scales))






  
     # titre_landraces_to_keep <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/landraces_to_keep.txt"
     # titre_landraces <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/landraces.txt"
     # titre_dissimilarities <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/landrace632_8741hap.dis"
     # dossier_sortie_liste <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/"
     # dossier_graphes <- "/home/adanguydesd/Documents/These_Alice/pipeline/160320/graphes/"
# titre_landraces_to_keep <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/landraces_PLINK_2.txt"

# titre_landraces <-  "/work/adanguy/these/pipeline/020820/sources/landraces.txt"      
# titre_dissimilarities <-"/work/adanguy/these/pipeline/amont/landrace632_8741hap.dis"     
# titre_landraces_to_keep <- "/work/adanguy/these/pipeline/020820/PLINK/landraces_PLINK_2.txt"
# dossier_sortie_liste <- "/work/adanguy/these/pipeline/020820/PLINK/"                     
# dossier_graphes <-"/work/adanguy/these/pipeline/020820/graphs/" 
 
# Thesholds
coeff_str <- 0.5 # minimal admixture coefficient to keep a landraces for further analysis, according K=4 STRUCTURE analysis
FDR <- 1e-3 # Treshold to identify higly related pairs of landraces



##################################################### Fonctions


#### Function 1 
# We have a matrix of dissimilarities = a N*N matrix, with dissimilarity measure for each pair of landraces
# Dissimilarity  ~ 1 - % of shared haplotypic alleles between a pair of landraces. Higly related landraces show dissimilarities ~ 0
# Function used to extract dissimilarities of landraces belonging to one population
dissimilarities_filtre <- function(dissimilarities, LINES_to_keep, LINES_to_delete){
  
  
  dissimilarities <- dissimilarities[ which(row.names(dissimilarities) %in% LINES_to_keep), which(colnames(dissimilarities) %in% LINES_to_keep)]
  
  
  
  # Exclure
  if(length(LINES_to_delete) >0){
    
    dissimilarities <- dissimilarities[- which(row.names(dissimilarities) %in% LINES_to_delete), -which(colnames(dissimilarities) %in% LINES_to_delete)]
    
  }
  
  # Output : filtered matrix
  return(dissimilarities)
  
}


##### Function 2  
# Convert dissimilarity data in a special format "adjency matrix" required to use the igraph package
# adjency matrix : = 1 when pairs of individual are higly related (dissimilarity above a treshold); = 0 otherwise
# Pairs of landraces exbiting 1 would be connected by a strait line in the graphical representation
convert_to_adjacency_matrix <- function(dissimilarities, seuil){
  
  # Chnage format
  dissimilarities <- as.matrix(dissimilarities)
  # Trick : dissimilarity > threshold are replaced by 100
  dissimilarities[dissimilarities > seuil] <- 100
  #  dissimality <= threshold are replace by 1
  dissimilarities[dissimilarities <= seuil] <- 1
  # Trick : dissimilarity > threshold are replaced by 0
  dissimilarities[dissimilarities == 100] <- 0
  
  # Sortie : an ajdency matrix (0 and 1)
  return(dissimilarities)
  
}


##### Function 3
# Produce a dataframe given connections between higly related landraces
# output : one row = a pair of connected landraces
obtain_connected_components <- function(adjency_matrix){
  
  # Function from igraph package
  connected_components_graph <- graph_from_adjacency_matrix(adjency_matrix, mode = "lower", weighted = T, diag = F, add.rownames = NA)
  # Extract connection (when adjency matrix = 1)
  connected_components <- as.data.frame(get.edgelist(connected_components_graph, names=T))
  
  # Sortie : data.frame avec un nb de lignes = toutes les paires d'individus apparentés, et deux colonnes contenant des ID
  return(connected_components)
  
}


##### Function 4 
# Count connection per landraces and sort landraces according to their number of connections
count_connections <- function(connected_components){
  
  # Names of connected individuals
  connections <- c(as.character(connected_components$V1), as.character(connected_components$V2))
  
  # If there is (still) connected l=individuals
  if(length(connections)>0){
    
    # Compt occurence of each name
    connections <- as.data.frame(table(connections))
    colnames(connections) <- c("ERGE","nb_connections")
    connections$ERGE <- as.character(connections$ERGE)
    # If two or more landraces exhibits the higher number of connections, if those landraces are not mixed, their is a risk to supress the landrace with the first ERGE according alpha-numeric sorting
    connections <- connections[sample(1:nrow(connections), size=nrow(connections)),]
    # Arrange the dataframe with higher connected landraces in the head
    connections <- connections[order(connections$nb_connections, decreasing=T),]
    
  }
  
  return(connections)
  
}

#### Function 5
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}


##### Step 1 : supress landraces with  missing data or hetozygotie
cat("\n\n input 1 : landraces with less than 5% heterozygotie and less than 10% missing data \n\n")
to_keep <- read.table(titre_landraces_to_keep, header=F, dec=".", sep="\t") 


head(to_keep)


to_keep <- to_keep %>%
  unlist() %>%
  as.vector()

length(to_keep)


##### Step 2 : removed hilgy admixed landraces

cat("\n\n input 2 : landraces passeport data \n\n")
l <- read.table(titre_landraces, header=T, dec=".", sep="\t") 

head(l)


l <- l%>%
  dplyr::select(-one_of("group","population","kept","popcolor")) %>%
  mutate(kept=ifelse(STRmemb4_max < coeff_str | ! LINE %in% to_keep | LINE =="KAMMAH", FALSE, TRUE)) 
# Kammah accession is supressed because it is duplicated in the data



table(l$kept)


##### Step 3 : Import dissimilarity matrix (~ distance matrix)
cat("\n\n input 3 : dissimilarity matrix of landraces \n\n")

nb_landraces=scan(titre_dissimilarities,skip=1,nlines=1, quiet = T)
dis.ssr=read.table(titre_dissimilarities,skip=3,fill=T,col.names=1:nb_landraces,nrows=(nb_landraces-1))
dis.ssr[1:5,1:5]
dis.ssr=dis.ssr[,2:ncol(dis.ssr)]
# Make a squared matrix N*N
dis.ssr=rbind.data.frame(NA,dis.ssr)
dis.ssr=cbind.data.frame(dis.ssr,NA)
# Landraces are named after their ERGE ID
colnames(dis.ssr) = row.names(dis.ssr) = l %>%
  arrange(Unit) %>%
  dplyr::select(ERGE) %>%
  unlist() %>%
  as.vector() %>%
  as.character()

dis.ssr[1:5,1:5]






##### Step 4 : hierarchical clustering

# Extract filtered landraces
to_keep <- l %>% filter(kept==T) %>%
  dplyr::select(ERGE) %>%
  unlist() %>%
  as.vector()

dist <- dis.ssr[which(row.names(dis.ssr) %in% to_keep ), which(colnames(dis.ssr) %in% to_keep)] %>%
  as.dist(dis.ssr)

# hierarchical clustering
arbre <- hclust(dist, method = "ward.D2")

K4 <- cutree(arbre, k=4) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(cluster=".", ERGE=rowname) %>%
  mutate(ERGE = as.character(ERGE))

K2 <- cutree(arbre, k=2) %>%
  as.data.frame() %>%
  rownames_to_column() %>% rename(group=".", ERGE=rowname) %>%
  mutate(ERGE = as.character(ERGE))

# Data formating
l <- l %>%
  mutate(ERGE = as.character(ERGE)) %>%
  full_join(K4, by="ERGE") %>%
  full_join(K2, by="ERGE") %>%
  mutate(group=paste0("G",group)) %>%
  mutate(cluster=paste0(group,"P",cluster)) %>%
  mutate(cluster=ifelse(cluster=="G2P3", "G2P1", cluster))%>%
  mutate(cluster=ifelse(cluster=="G2P4", "G2P2", cluster))%>%
  mutate(cluster=ifelse(kept=="FALSE", NA, cluster)) %>%
  mutate(population=case_when(
    cluster=="G1P1"~ "WA",
    cluster=="G1P2"~"EA",
    cluster=="G2P1"~"WE",
    cluster=="G2P2"~"EE")) %>%
  mutate(group=case_when(
    group=="G1"~"ASIA",
    group=="G2"~"EUROPE")) %>%
  mutate(popcolor=case_when(
    cluster=="G1P1"~"#ABDDA4" ,
    cluster=="G1P2"~"#2B83BA",
    cluster=="G2P2"~"#FDAE61",
    cluster=="G2P1"~"#D7191C")) %>%
  mutate(STRcolor8=case_when(
    STRarea8=="NWE" ~ "#D7191C",
    STRarea8=="IBP" ~ "#FEE08B",
    STRarea8=="SEE" ~ "#FDAE61",
    STRarea8=="MED" ~ "#F46D43",
    STRarea8=="CAU" ~ "#E6F598",
    STRarea8=="CAA" ~ "#ABDDA4",
    STRarea8=="INP" ~ "#66C2A5",
    STRarea8=="SEA" ~ "#2B83BA"))  %>%
  dplyr::select(-one_of("cluster")) %>%
  arrange(LINE)

table(l$population, l$STRarea8)


##### Step 5 : supress higly related landraces


LINES_to_delete_final <- as.character()


# Within each population, the objectiv is to identify higly related landraces and supress them
for (c in unique(l$population) %>% na.omit() %>% as.vector()){
  
  
  
  # Names of landraces
  LINES_to_keep <- l %>% filter(population== !!c) %>%
    filter(kept == T) %>%
    dplyr::select("ERGE") %>%
    unlist() %>% as.vector()
  
  LINES_to_delete <- as.character()
  
  
  
  
  
  #### Identification of higly related landraces
  
  # Need of a vector, used to draw a distribution. Here it is a vector of dissimilarity, expected to be close of 0 when a pair of landraces are similar
  # One value = one pair of landraces
  vecteur_dissimilarities_cluster <- dis.ssr %>% dissimilarities_filtre(., LINES_to_keep = LINES_to_keep, LINES_to_delete = LINES_to_delete) %>%
    unlist() %>%
    as.vector() %>%
    na.omit()
  
  
  
  

  # Compute robust mean and standard deviation with robust linear modelisation
  # Iterativ procedur reducing progressively weights of outliers when computing mean and variance of a fitted Normal distribution
  # Package MASS
  modelisation <- rlm(vecteur_dissimilarities_cluster~1)
  moyenne <- as.vector(modelisation$coefficients)
  ecart_type <- as.vector(modelisation$s)
  
  # Compute for value the probability to belong to the Normale distributon
  # Lower tail =T used to target higly related pairs of landraces, expected to have low value of dissimilarity 
  proba_brutes <- pnorm(vecteur_dissimilarities_cluster, mean=moyenne, sd=ecart_type, lower.tail = T)
  
  # Correcting for multiple testing. Compute False Discovery Rate with qvalue package
  proba_corrigees <- qvalue(proba_brutes)$qvalues
  
  # Pairs of landraces exibiting a FDR <= 1e-3 are considered as higly related
  seuil <- max(vecteur_dissimilarities_cluster[which(proba_corrigees <= FDR)],0)
  
  
  # Supress higly related pairs of landraces progressively, starting with landraces exibitng the higher number of relationship, until no close relationship remains

    # (la variable sortie indique si le data.frame connected_comonents à l'étape contient encore des individus à supprimer
  # On s'arrête quand sortie = 1 car cela signfie qu'on a supprimé un des deux individus de la dernière paire apparentée
  # liste des individus à supprimer (vecteur))
  
  connections <-  dis.ssr %>% dissimilarities_filtre(., LINES_to_keep = LINES_to_keep, LINES_to_delete = LINES_to_delete) %>%
    convert_to_adjacency_matrix(dissimilarities=., seuil=seuil) %>%
    obtain_connected_components(.) %>%
    count_connections(.)
  
  
  while (length(connections)>0){
    
    
    # If theire is still higly related individuals, extract the one exibiting the more relationships
    LINE <- unique(connections$ERGE)[1]
    # Premier individu = celui avec le + de connections
    
    # Add its ID to the list of ID of landraces to supress
    LINES_to_delete <- c(LINES_to_delete, LINE)
    
    connections <-  dis.ssr %>% dissimilarities_filtre(., LINES_to_keep = LINES_to_keep, LINES_to_delete = LINES_to_delete) %>%
      convert_to_adjacency_matrix(dissimilarities=., seuil=seuil) %>%
      obtain_connected_components(.) %>%
      count_connections(.)
    
    
    
    
  }
  
  
  
  # Final list of landraces to supress
  LINES_to_delete_final <- c(LINES_to_delete_final, LINES_to_delete)
  
  
  
  
##### Step 6 : graphical representation of relatedness
  
  titre_graphe <- paste0(dossier_graphes,c,".png")
  png(titre_graphe)
  layout(matrix(c(1,3,2,3), 2, 2, byrow = TRUE),     widths=c(1,2), heights=c(1,1))
  par(mar = c(2,2,6,1), xpd=T)
  
  # Dissimilarities distribution
  hist(vecteur_dissimilarities_cluster, main = "Distance", xlab="", xlim=c(0,1), freq = F, breaks = 20)
  segments(x0=seuil,y0=0,x1=seuil,y1=6, col="red", lty=2)
  legend(x=0.7, y=6, c(round(seuil, digits=2)), col="red", lty=2, bty="n", cex=0.8, title="Treshold")
  xx=seq(0,1,0.01) # gamme de valeurs pour l'axe des abscisses
  d1=dnorm(xx, mean = moyenne, sd = ecart_type)
  lines(xx, d1, lwd=1, col="blue")
  
# Raw probability to belong to the Normal distribution
    hist(proba_brutes, xlab ="", xlim=c(0,1), main="Raw probabilities", breaks = 20)
  
  # Graoh of connected individuals
  vecteur <- dis.ssr %>% dissimilarities_filtre(., LINES_to_keep = LINES_to_keep, LINES_to_delete = as.character()) %>%
    row.names() %in% 
    LINES_to_delete 
  vecteur[vecteur==TRUE] <- "red"
  vecteur[vecteur==FALSE] <- "blue"
  adjency_matrix <- dis.ssr %>% dissimilarities_filtre(., LINES_to_keep = LINES_to_keep, LINES_to_delete = as.character()) %>%
    convert_to_adjacency_matrix(dissimilarities=., seuil=seuil)
  test <- graph_from_adjacency_matrix(adjency_matrix, mode = "lower", weighted = T, diag = F)
  plot(test,vertex.size=3, vertex.color=vecteur, main="Connected Components", vertex.label=NA)
  
  main=gsub("landraces_","",gsub(pattern=c(".txt"), replacement = "", c))
  title(main, line = -2, outer = TRUE, font.main= 4, cex.main=2, sub = paste0(length(LINES_to_keep), " landraces (",length(LINES_to_delete)," discarded)"))
  legend(0.7,1, c("discarded"), col="red", pch=19, bg="transparent")
  
  
  box(which = "outer")
  
  
  dev.off()
  
}


#### Step 7 : data output

l<- l %>% mutate(kept = ifelse(ERGE%in% LINES_to_delete_final | kept ==F, "FALSE","TRUE")) %>%
  arrange(LINE) %>%
  mutate(popcolor=ifelse(kept==F,"black",popcolor))



as.data.frame(table(l$population, l$STRgp4)) %>% 
  mutate(ntot=sum(Freq)) %>%
  arrange(desc(Freq)) %>%
  group_by(Var1) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(n=sum(Freq)) %>%
  summarise(prop_consistentK4=round(unique(n/ntot),2))




as.data.frame(table(l$population, l$STRgp8)) %>% 
  mutate(ntot=sum(Freq)) %>%
  arrange(desc(Freq)) %>%
  group_by(Var2) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(n=sum(Freq)) %>%
  summarise(prop_consistentK8=round(unique(n/ntot),2))


cat("\n\n Output 1 : passeport data for landraces + hierarchical clustering and filtering on relatedness \n\n")
write.table(l, titre_landraces, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
head(l)

for (c in unique(l$population) %>% na.omit() %>% as.vector()){
  
  
  titre=paste0(dossier_sortie_liste, c,".txt")
  
  l2 <- l %>% filter(population==!!c & kept==T) %>%
    mutate(V1=1) %>%
    dplyr::select(V1, LINE) %>%
    arrange(LINE)
  
  cat("\n\n Output 2 : final list of landraces per population WE, EE, WA, EA \n\n")
  write.table(l2, titre, col.names = F, row.names = F, dec=".", sep="\t", quote=F)
  head(l2)
}





# 
# #### Step 8 : Graphical representation of hierarchical clustering + distance matrix + STRUCTURE K=4 and K=8
# 
# 
# arbre2 <- as.dendrogram(arbre) 
# ordre <- labels(arbre2)
# 
# 
# 
# l2 <- l %>% filter(!is.na(population))
# l2 <- l2[match(ordre,l2$ERGE),]
# 
# l3 <- l2 %>% filter(kept==T)
# 
# 
# 
# arbre2 <- arbre2 %>% color_branches(., k=4, col = unique(l3$popcolor)) %>%
#   set("leaves_pch", 19) %>%
#   set("leaves_col",l2$popcolor) %>%
#   set("leaves_cex", 0.5) %>%
#   set("branches_lwd", 1)
# 
# dis.ssr2 <- makeSymm(dis.ssr)
# diag(dis.ssr2) <- 0
# dis.ssr2 <- dis.ssr2[match(l2$ERGE,colnames(dis.ssr)),match(l2$ERGE,colnames(dis.ssr)) ]
# 
# 
# d <- dis.ssr2 %>%
#   rownames_to_column('ERGE1') %>%
#   gather(ERGE2, Distance, -ERGE1) %>%
#   mutate(ERGE1 = factor(ERGE1, levels=unique(ERGE1))) %>%
#   mutate(ERGE2 = factor(ERGE2, levels=unique(ERGE2))) %>%
#   mutate(pos1=as.factor(as.numeric(ERGE1))) %>%
#   mutate(pos2=rev(as.factor(as.numeric(ERGE2))))
# 
# 
# 
# # Admixture matrix for each landrace (K=4)
# qSTR4 <- l2 %>% dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
#   dplyr::select((-one_of("STRmemb4_max"))) %>%
#   gather(memb, value, -ERGE) %>%
#   mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
#   mutate(pos=as.factor(as.numeric(ERGE))) %>%
#   mutate(memb=str_replace_all(memb,"STRmemb4_",""))
# 
# 
# t4 <- l3 %>% dplyr::select(popcolor, starts_with("STRmemb4_")) %>%
#   dplyr::select((-one_of("STRmemb4_max"))) 
# 
# t4$memb <- paste0("STRmemb4_",as.vector(apply(t4[,-c(1)], 1, function(x) which.max(x)))) 
# 
# 
# t4 <- as.data.frame(table(t4$popcolor, t4$memb)) %>%
#   arrange(desc(Freq )) %>%
#   slice(1:4) %>%
#   rename(popcolor=Var1, memb=Var2) %>%
#   arrange(memb) %>%
#   mutate(memb=str_replace_all(memb,"STRmemb4_",""))
# 
# 
# 
# qSTR4 <- qSTR4 %>%
#   inner_join(t4, by="memb") 
# 
# # head(qSTR4)
# 
# 
# 
# 
# 
# 
# # Admixture matrix for each landrace (K=8)
# qSTR8 <- l2 %>% dplyr::select(ERGE, starts_with("STRmemb8_")) %>%
#   dplyr::select((-one_of("STRmemb8_max"))) %>%
#   gather(memb, value, -ERGE) %>%
#   mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
#   mutate(pos=as.factor(as.numeric(ERGE)))
# 
# t8 <- l3 %>% dplyr::select(STRcolor8, STRarea8, starts_with("STRmemb8_")) %>%
#   dplyr::select((-one_of("STRmemb8_max"))) 
# 
# t8$memb <- paste0("STRmemb8_",as.vector(apply(t8[,-c(1,2)], 1, function(x) which.max(x)))) 
# 
# 
# t8 <- as.data.frame(table(t8$STRcolor8, t8$STRarea8, t8$memb)) %>%
#   arrange(desc(Freq )) %>%
#   slice(1:8) %>%
#   rename(STRcolor8=Var1, STRarea8=Var2, memb=Var3) %>%
#   arrange(memb) 
# 
# 
# qSTR8 <- qSTR8 %>%
#   inner_join(t8, by="memb") %>%
#   mutate(memb=str_replace_all(memb,"STRmemb8","K8")) %>%
#   mutate(memb=as.factor(memb)) %>%
#   mutate(STRarea8=factor(STRarea8, levels=t8$STRarea8))
# 
# 
# # head(qSTR8)
# # head(t8)
# 
# 
# titre_abre <- paste0(dossier_graphes,"tree.png")
# png(titre_abre)
# plot(arbre2, main="", leaflab="none", axes=FALSE, ylab="")
# rect.dendrogram(arbre2, k=8, border="grey", lwd=0.5, xpd=F)
# legend(547,3,c(rev(unique(l3$population)),"discarded"), y.intersp=4, col=c(rev(unique(l3$popcolor)), "black"),title = "Population", text.font=2, pch=15, bty = "n", cex=0.4, pt.cex =0.9)
# dev.off()
# 
# par(mar=c(0,0,0,1.5), oma=c(0,0,0,0), xpd=TRUE, font=2, family="serif")
# plot(arbre2, main="", leaflab="none", axes=FALSE, ylab="")
# legend(547,3,c(rev(unique(l3$population)),"discarded"), y.intersp=4, col=c(rev(unique(l3$popcolor)), "black"),title = "Population", text.font=2, pch=15, bty = "n", cex=0.4, pt.cex =0.9)
# 
# p1 <- recordPlot()
# 
# 
# 
# p2 <- ggplot(d, aes(pos1, pos2)) +
#   geom_tile(aes(fill = Distance)) + 
#   scale_fill_gradient(low = "white", high = "red") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   theme(plot.margin = unit(c(0,0.1,0,0.5), "cm")) +
#   theme(legend.text = element_text(size = 5)) +
#   theme(legend.key.size = unit(0.5,"line")) +
#   theme(text=element_text(family="serif", face="bold")) +
#   theme(legend.title=element_text(size=5)) +
#   theme(legend.box.spacing=unit(0.2, 'cm'))
# 
# 
# p3 <- ggplot(data=qSTR4, aes(x=pos, y=value, fill=memb, col=memb) ) + 
#   geom_bar(stat="identity") +
#   theme_void ()+
#   scale_colour_manual("K=4", values=as.character(t4$popcolor)) +
#   scale_fill_manual("K=4",  values=as.character(t4$popcolor)) +
#   theme(plot.margin = unit(c(0,0.5,0,0.6), "cm")) +
#   theme(legend.text = element_text(size = 5)) +
#   theme(legend.key.size = unit(0.5,"line")) +
#   theme(text=element_text(family="serif", face="bold")) +
#   theme(legend.title=element_text(size=5)) +
#   scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
#   theme(legend.box.spacing=unit(0.4, 'cm'))
# 
# 
# 
# 
# 
# 
# p4 <- ggplot(data=qSTR8, aes(x=pos, y=value, col=STRarea8, fill=STRarea8) ) + 
#   geom_bar(stat="identity") +
#   theme_void ()+
#   scale_colour_manual("K=8", values=as.character(t8$STRcolor8)) +
#   scale_fill_manual("K=8", values=as.character(t8$STRcolor8)) +
#   theme(plot.margin = unit(c(0,0.15,0,0.6), "cm")) +
#   theme(legend.text = element_text(size = 5)) +
#   theme(legend.key.size = unit(0.5,"line")) +
#   theme(text=element_text(family="serif", face="bold")) +
#   theme(legend.title=element_text(size=5)) +
#   scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
#   theme(legend.box.spacing=unit(0.4, 'cm'))
# 
# ggarrange(p1,p2,p3,p4, ncol = 1, nrow = 4, labels = c("A","B","C","D"))
# 
# 
# titre_figure <- paste0(dossier_graphes,"landraces.png")
# 
# figure <- ggarrange(p1,p2,p3,p4, ncol = 1, nrow = 4, labels = c("A","B","C","D"))
# 
# cat("\n\n Graph 1 : hierarchical clusering + distance matrix + STRUCTURE K=4 and K=8 \n\n")
# ggsave(titre_figure, plot = last_plot(), width = 7.29, height = 4.5, units = c("in"))
# 



sessionInfo()