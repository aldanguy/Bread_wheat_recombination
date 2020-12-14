  # Drawn Figure 1 : FST tree, admixture graphes for K=4 and K=8 and FST matrix
  
  
  
  Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(qgraph))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gplots))
#suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gridGraphics))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(viridis))



cat("\n landraces_graphs.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_landraces <- variables[1]
titre_tree <- variables[2]
titre_matrice <- variables[3]
titre_dissimilarities <- variables[4]
titre_graphe1 <-variables[5]
titre_graphe2 <-variables[6]
titre_matrice_FST_png<-variables[7]
titre_hclust_png <- variables[8]
titre_tabs2 <- variables[9]


    titre_landraces <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/landraces.txt"
    titre_tree <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/hapflk_tree_genome.txt"
    titre_matrice <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/pairwise_FST_matrix_haplotypic_blocks.txt"
    titre_graphe2 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure1.tiff"
    titre_matrice_FST_png <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/landraces_FST_matrix_haplotypes.png"
    titre_hclust_png <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/hclust4_png.png"
    titre_graphe1 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/landraces.png"
   titre_dissimilarities <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/amont/landrace632_8741hap.dis"
   


# titre_landraces <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/landraces.txt"
# titre_matrice <-  "/home/adanguydesd/Documents/These_Alice/pipeline/020820/pairwise_FST_matrix_haplotypic_blocks.txt"
# titre_graphe <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/graphes/landraces_2.png"
# titre_tree <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/hapflk_tree_genome.txt"
# titre_dissimilarities <- "/home/adanguydesd/Documents/These_Alice/pipeline/amont/landrace632_8741hap.dis"
# titre_hclust_png <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/graphes/landraces_hclust4.png"
# titre_matrice_FST_png <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/graphes/landraces_FST_matrix_haplotypes.png"
# titre_graphe1 <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/graphes/landraces1.png"
# titre_graphe2 <- "/home/adanguydesd/Documents/These_Alice/pipeline/020820/graphes/landraces2.png"


makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}


color4 <- data.frame(pop=c("NWE","SEE","CAA","SEA"), col=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA"), gp=c("3","4","1","2"))
color8 <- data.frame(pop=c("NWE","MED","SEE","IBP","CAU","CAA","INP","SEA"), col=c("#D7191C","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#2B83BA"), gp=c("5","6","8","7","2","1","3","4"))



cat("\n\n Input 1 : landraces data \n")
l <- fread(titre_landraces) 
head(l)



cat("\n\n Input 2 : tree of genome wide Fst from Reynolds distance \n")
tree <- read.tree(titre_tree)
tree




cat("\n\n Input 3 : Fst matrix \n")
matrice <- read.table(titre_matrice, row.names = 1)
matrice <- as.matrix(matrice)
matrice

# matrice <- matrice %>% column_to_rownames("V1")
# colnames(matrice) <- row.names(matrice)
# col.ordre <- c("WE","EE","WA","EA")
# matrice <- matrice[col.ordre,col.ordre]


cat("\n\n Input 4 : distance matrix \n\n")
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
dis.ssr2 <- makeSymm(dis.ssr)
diag(dis.ssr2) <- 0

##### 1) hclust tree + distance matrix + K=4 & K=8

to_keep <- l %>% filter(!is.na(population)) %>%
  dplyr::select(ERGE) %>%
  unlist() %>%
  as.vector()

dist <- dis.ssr[which(row.names(dis.ssr) %in% to_keep ), which(colnames(dis.ssr) %in% to_keep)] %>%
  as.dist(dis.ssr)



arbre <- hclust(dist, method = "ward.D2")
arbre2 <- as.dendrogram(arbre) 


arbre2 <- all_couple_rotations_at_k(arbre2, k = 3)[[2]]
arbre2 <- all_couple_rotations_at_k(arbre2, k = 4)[[2]]
arbre2 <- all_couple_rotations_at_k(arbre2, k = 5)[[2]]
arbre2 <- all_couple_rotations_at_k(arbre2, k = 7)[[2]]

ordre <- labels(arbre2)


l2 <- l %>% filter(!is.na(population))
l2 <- l2[match(ordre,l2$ERGE),]
l3 <- l2 %>% filter(kept==T)


dis.ssr2 <- dis.ssr2[match(l2$ERGE,colnames(dis.ssr)),match(l2$ERGE,colnames(dis.ssr)) ]
d <- dis.ssr2 %>%
  rownames_to_column('ERGE1') %>%
  gather(ERGE2, Distance, -ERGE1) %>%
  mutate(ERGE1 = factor(ERGE1, levels=unique(ERGE1))) %>%
  mutate(ERGE2 = factor(ERGE2, levels=unique(ERGE2))) %>%
  mutate(pos1=as.factor(as.numeric(ERGE1))) %>%
  mutate(pos2=rev(as.factor(as.numeric(ERGE2))))


qSTR4 <- l2 %>% dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>%
  gather(memb, value, -ERGE) %>%
  mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
  mutate(pos=as.factor(as.numeric(ERGE))) %>%
  mutate(memb=str_replace_all(memb,"STRmemb4_",""))

t4 <- l3 %>% dplyr::select(popcolor, starts_with("STRmemb4_")) %>%
  dplyr::select((-one_of("STRmemb4_max"))) 
t4$memb <- paste0("STRmemb4_",as.vector(apply(t4[,-c(1)], 1, function(x) which.max(x)))) 
t4 <- as.data.frame(table(t4$popcolor, t4$memb)) %>%
  arrange(desc(Freq )) %>%
  slice(1:4) %>%
  rename(popcolor=Var1, memb=Var2) %>%
  arrange(memb) %>%
  mutate(memb=str_replace_all(memb,"STRmemb4_",""))



qSTR4 <- qSTR4 %>%
  inner_join(t4, by="memb") 

# Admixture matrix for each landrace (K=8)
qSTR8 <- l2 %>% dplyr::select(ERGE, starts_with("STRmemb8_")) %>%
  dplyr::select((-one_of("STRmemb8_max"))) %>%
  gather(memb, value, -ERGE) %>%
  mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
  mutate(pos=as.factor(as.numeric(ERGE)))

t8 <- l3 %>% dplyr::select(STRcolor8, STRarea8, starts_with("STRmemb8_")) %>%
  dplyr::select((-one_of("STRmemb8_max"))) 
t8$memb <- paste0("STRmemb8_",as.vector(apply(t8[,-c(1,2)], 1, function(x) which.max(x)))) 
t8 <- as.data.frame(table(t8$STRcolor8, t8$STRarea8, t8$memb)) %>%
  arrange(desc(Freq )) %>%
  slice(1:8) %>%
  rename(STRcolor8=Var1, STRarea8=Var2, memb=Var3) %>%
  arrange(memb) 

qSTR8 <- qSTR8 %>%
  inner_join(t8, by="memb") %>%
  mutate(memb=str_replace_all(memb,"STRmemb8","K8")) %>%
  mutate(memb=as.factor(memb)) %>%
  mutate(STRarea8=factor(STRarea8, levels=t8$STRarea8))


arbre2 <- arbre2 %>% color_branches(., k=4, col = c("#D7191C", "#FDAE61","#ABDDA4","#2B83BA")) %>%
  set("leaves_pch", 19) %>%
  set("leaves_col",l2$popcolor) %>%
  set("leaves_cex", 2) %>%
  set("branches_lwd", 5)


png(titre_hclust_png, width = 50, height = 10, units = 'in',res = 300)
par(mai=c(0,0,0,2), xpd=TRUE, font=2, family="serif")
plot(arbre2, main="", leaflab="none", axes=FALSE, ylab="")
legend(540,4,c("WE","EE","WA","EA","discarded"), inset=c(10,0), col=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA", "black"),title = "Population", text.font=2, pch=15, bty = "n", cex=3, pt.cex =5)

dev.off()


p2 <- ggplot(d, aes(pos1, pos2)) +
  geom_tile(aes(fill = Distance)) + 
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin = unit(c(0,0.1,0,0.5), "cm")) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.5,"line")) +
  theme(text=element_text(family="serif", face="bold")) +
  theme(legend.title=element_text(size=5)) +
  theme(legend.box.spacing=unit(0.2, 'cm'))


p3 <-   qSTR4 %>%
  mutate(memb=factor(memb, levels = c("3","4","1","2"))) %>%
  ggplot(aes(x=pos, y=value, fill=memb, col=memb) ) + 
  geom_bar(stat="identity") +
  theme_void ()+
  scale_colour_manual("K=4", values=as.character(color4$col), labels=color4$pop) +
  scale_fill_manual("K=4",  values=as.character(color4$col), labels=color4$pop) +
  theme(plot.margin = unit(c(0,0.2,0,0.6), "cm")) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.5,"line")) +
  theme(text=element_text(family="serif", face="bold")) +
  theme(legend.title=element_text(size=5)) +
  scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
  theme(legend.box.spacing=unit(0.4, 'cm'))



p4 <- qSTR8 %>%
  mutate(STRarea8=factor(STRarea8, levels = color8$pop)) %>%
  ggplot(aes(x=pos, y=value, col=STRarea8, fill=STRarea8) ) + 
  geom_bar(stat="identity") +
  theme_void ()+
  scale_colour_manual("K=8", values=as.character(color8$col), labels=color8$pop) +
  scale_fill_manual("K=8", values=as.character(color8$col), labels=color8$pop) +
  theme(plot.margin = unit(c(0,0.15,0,0.6), "cm")) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.5,"line")) +
  theme(text=element_text(family="serif", face="bold")) +
  theme(legend.title=element_text(size=5)) +
  scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
  theme(legend.box.spacing=unit(0.4, 'cm'))
g <- ggdraw() +
   draw_plot(p2, x=0, y=0.66,1, 0.33) +
   draw_plot(p3,x=0,y=0.33,1,height=0.33) +
   draw_plot(p4,x=0,y=0,1,height=0.33) +
  draw_plot_label(c("B", "C", "D"), c(0.93,0.93,0.93), c(1,0.6,0.33), size = 15)

graph1 <- ggdraw() +
  draw_image(titre_hclust_png, x=0, y=0.73, width=1, height=0.25, scale=1, clip = T) +
  draw_plot(g,x=0,y=0,1,0.75)+
  draw_plot_label(c("A", ""), c(0.93,1), c(1,0), size = 15)

cat("\n\n Output 1 : graph including hclust, distance matrice K=4+8")
graph1
ggsave(titre_graphe1, graph1)




#### 2) FST tree + K=4+8


l <- l %>%
  filter(kept==T) 

t4 <- l %>% dplyr::select(popcolor, starts_with("STRmemb4_")) %>%
  dplyr::select((-one_of("STRmemb4_max"))) 
t4$memb <- paste0("STRmemb4_",as.vector(apply(t4[,-c(1)], 1, function(x) which.max(x)))) 
t4 <- as.data.frame(table(t4$popcolor, t4$memb)) %>%
  arrange(desc(Freq )) %>%
  slice(1:4) %>%
  rename(popcolor=Var1, memb=Var2) %>%
  arrange(memb) %>%
  mutate(memb=str_replace_all(memb,"STRmemb4_",""))
# 
# t8 <- l %>% dplyr::select(STRcolor8, starts_with("STRmemb8_")) %>%
#   dplyr::select((-one_of("STRmemb8_max"))) 
# t8$memb <- paste0("STRmemb8_",as.vector(apply(t8[,-c(1)], 1, function(x) which.max(x)))) 
# t8 <- as.data.frame(table(t8$STRcolor8, t8$memb)) %>%
#   arrange(desc(Freq )) %>%
#   slice(1:8) %>%
#   rename(popcolor=Var1, memb=Var2) %>%
#   arrange(memb) %>%
#   mutate(memb=str_replace_all(memb,"STRmemb8_",""))



color4 <- c("#D7191C","#FDAE61","#ABDDA4","#2B83BA")

#ABDDA4 = green
#D7191C = blue
#FDAE61 = red
#2B83BA  = orange

## WE




WE4 <- l %>% filter(population=="WE")%>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  mutate(value_dom=STRmemb4_4+STRmemb4_1+STRmemb4_2) %>%
  mutate(value_dom2=STRmemb4_4) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE, -value_dom, -value_dom2) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>%
  mutate(memb=case_when(memb==3~"NWE",
                        memb==4~"SEE",
                        memb==1~"CAA",
                        memb==2~"SEA")) %>%
  mutate(memb=factor(memb, levels=c("NWE","SEE","CAA","SEA"))) %>%
  arrange(value_dom, desc(value_dom2)) %>%
  mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
  mutate(pos=as.factor(as.numeric(ERGE))) %>%
  mutate(memb=factor(memb, levels=c("NWE","SEE","CAA","SEA"))) 


WEg4 <- WE4  %>%
  ggplot(aes(x=pos, y=value, fill=memb, col=memb) ) + 
  geom_bar(stat="identity", position="fill") +
  theme_void ()+
  scale_colour_manual("K=4", values=color4) +
  scale_fill_manual("K=4",  values=color4) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.5,"line")) +
  theme(text=element_text(family="serif", face="bold")) +
  theme(legend.title=element_text(size=5)) +
  scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
  theme(legend.box.spacing=unit(1, 'cm')) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.text=element_text(size=10))+
  theme(plot.margin = unit(c(1,0.05,0,0.05), "cm"))

WEg4






##### EE



EE4 <- l %>% filter(population=="EE")%>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  mutate(value_dom=STRmemb4_3) %>%
  mutate(value_dom2=STRmemb4_1) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE, -value_dom, -value_dom2) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>%
  mutate(memb=case_when(memb==3~"NWE",
                        memb==4~"SEE",
                        memb==1~"CAA",
                        memb==2~"SEA")) %>%
  mutate(memb=factor(memb, levels=c("SEE","NWE","CAA","SEA"))) %>%
  arrange(desc(value_dom), value_dom2) %>%
  mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
  mutate(pos=as.factor(as.numeric(ERGE))) %>%
  mutate(memb=factor(memb, levels=c("NWE","SEE","CAA","SEA"))) 




EEg4 <- EE4 %>%
  ggplot(aes(x=pos, y=value, fill=memb, col=memb) ) + 
  geom_bar(stat="identity", position="fill") +
  theme_void ()+
  scale_colour_manual("K=4", values=color4) +
  scale_fill_manual("K=4",  values=color4) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.5,"line")) +
  theme(text=element_text(family="serif", face="bold")) +
  theme(legend.title=element_text(size=5)) +
  scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
  theme(legend.box.spacing=unit(1, 'cm')) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.text=element_text(size=10))+
  theme(plot.margin = unit(c(1,0.05,0,0.05), "cm"))


EEg4


#### WA



WA4 <- l %>% filter(population=="WA")%>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  mutate(value_dom=STRmemb4_4) %>%
  mutate(value_dom2=STRmemb4_2) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE, -value_dom, -value_dom2) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>%
  mutate(memb=case_when(memb==3~"NWE",
                        memb==4~"SEE",
                        memb==1~"CAA",
                        memb==2~"SEA")) %>%
  mutate(memb=factor(memb, levels=c("SEE","NWE","CAA","SEA"))) %>%
  arrange(desc(value_dom), value_dom2) %>%
  mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
  mutate(pos=as.factor(as.numeric(ERGE))) %>%
  mutate(memb=factor(memb, levels=c("NWE","SEE","CAA","SEA"))) 



WAg4 <- WA4 %>%
  ggplot(aes(x=pos, y=value, fill=memb, col=memb) ) + 
  geom_bar(stat="identity", position="fill") +
  theme_void ()+
  scale_colour_manual("K=4", values=color4) +
  scale_fill_manual("K=4",  values=color4) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.5,"line")) +
  theme(text=element_text(family="serif", face="bold")) +
  theme(legend.title=element_text(size=5)) +
  scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
  theme(legend.box.spacing=unit(1, 'cm')) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.text=element_text(size=10))+
  theme(plot.margin = unit(c(1,0.05,0,0.05), "cm"))

WAg4

### EA



EA4 <- l %>% filter(population=="EA")%>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  mutate(value_dom=STRmemb4_4+STRmemb4_3+STRmemb4_1) %>%
  mutate(value_dom2=STRmemb4_1) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE, -value_dom, -value_dom2) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>%
  mutate(memb=case_when(memb==3~"NWE",
                        memb==4~"SEE",
                        memb==1~"CAA",
                        memb==2~"SEA")) %>%
  mutate(memb=factor(memb, levels=c("SEE","NWE","CAA","SEA"))) %>%
  arrange(desc(value_dom), desc(value_dom2)) %>%
  mutate(ERGE = factor(ERGE, levels=unique(ERGE))) %>%
  mutate(pos=as.factor(as.numeric(ERGE))) %>%
  mutate(memb=factor(memb, levels=c("NWE","SEE","CAA","SEA"))) 


EAg4 <- EA4 %>%
  ggplot(aes(x=pos, y=value, fill=memb, col=memb) ) + 
  geom_bar(stat="identity", position="fill") +
  theme_void ()+
  scale_colour_manual("K=4", values=color4) +
  scale_fill_manual("K=4",  values=color4) +
  theme(legend.text = element_text(size = 5)) +
  theme(legend.key.size = unit(0.5,"line")) +
  theme(text=element_text(family="serif", face="bold")) +
  theme(legend.title=element_text(size=5)) +
  scale_y_continuous(limits=c(0,0.95),oob = rescale_none) +
  theme(legend.box.spacing=unit(1, 'cm')) +
  theme(legend.title=element_text(size=15)) +
  theme(legend.text=element_text(size=10))+
  theme(plot.margin = unit(c(1,0.05,0,0.05), "cm"))


EAg4




K4 <- ggarrange(WEg4, EEg4, WAg4, EAg4, ncol=4, common.legend = T, legend = "none")


# K <- ggarrange(K4, nrow=1, labels = c("Population structure"), font.label = list(size=20))
K <- ggarrange(K4, nrow=1)










#### FST tree







# WEEE=0.01513686
# WEWA=0.03652387
# WEEA=0.08456574
# EEWA=0.03345497
# EEEA=0.08445896
# WAEA=0.04125292
# 
# matrice <- matrix(data=c(0,WEEE,WEWA,WEEA,
#                          WEEE,0,EEWA,EEEA,
#                          WEWA,EEWA,0,WAEA,
#                          WEEA,EEEA,WAEA,0), nrow=4, byrow=T)
# row.names(matrice) <- c("WE","EE","WA","EA")
# colnames(matrice) <- row.names(matrice)



# test <- nj(matrice)
#test <- as.phylo(hclust(dist(matrice)))
test <- tree
d = data.frame(node=1:7, color=c("#ABDDA4","#2B83BA", "#FDAE61","#D7191C","black","black","black"))
t <- ggtree(test, size=2) %<+% d + aes(color=I(color)) +
  theme_tree() +
  theme(legend.position="none") +
  scale_y_reverse() +
  scale_x_reverse() +
  coord_flip() +
  geom_tiplab(geom="label", offset=-0.005, size=10, align=F, hjust = 0.5, fontface="bold") +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) 



t <- ggtree::rotate(t, 6) 

# geom_treescale(x=0, y=-4, offset = 0.1) +




arbre <- ggdraw() +
  draw_plot(t, x=0, .48, 0.9, .5) +
  draw_plot(K, 0, 0, 1, .5) +
  draw_plot_label(c("Population tree", "Population structure"), c(-0.09, -0.12), c(1, 0.5), size = 20)




### FST matrix
n=4
t =10
  # 
dev.off()
png(titre_matrice_FST_png, units="in", width = 24, height = 24, res=200)

par(xpd = TRUE)
corrplot(round(matrice*100,digits=2), method="color", col=c(rep("gray94", times=11000),
                                                            rep(viridis(t)[n+6], times=2000),
                                                            rep(viridis(t)[n+5], times=1100),
                                                            rep(viridis(t)[n+4], times=700),
                                                            rep(viridis(t)[n+3], times=5200),
                                                            rep(viridis(t)[n+2], times=100),
                                                            rep(viridis(t)[n+1], times=1)
                                                            
),  
type="upper", 
addCoef.col = "black", # Ajout du coefficient de corrélation
tl.col=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA"), tl.srt=45, sig.level = 0.01, insig = "blank", 
# Cacher les coefficients de corrélation sur la diagonale
diag=T , is.corr = F, shade.lwd=0,number.digits=1, cl.pos = "n", tl.cex=10, number.cex=10,mar = c(0, 0, 0, 0), font.label = list(style = "bold"))
dev.off()
p <- readPNG(titre_matrice_FST_png)


graphe <- ggdraw() +
  draw_image(titre_matrice_FST_png, 0.7, 0.2, 0.3, 1) +
  draw_plot(arbre, x=0, 0, 0.7, 1) +
  draw_plot_label(c("Fst matrix (%)", ""), c(0.7, 0.7), c(1, 1), size = 20)

graphe
cat("\n\n Output 1 : landraces Fst and admixture \n")
tiff(titre_graphe2, units="in", width = 15, height = 9, res=200, compression = "lzw")
graphe
dev.off()
#ggsave(titre_graphe2, graphe, width = 32, height=20, units="cm")



#ABDDA4 = green
#D7191C = blue
#FDAE61 = red
#2B83BA = orange

l2 <- l %>% 
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>% 
  group_by(population, country) %>%
  summarise(n=n()) %>% 
  arrange(population, desc(n), country) %>% 
  group_by(population) %>% 
  mutate(total=sum(n)) %>%
  rename(group=population)




write.table(l2, titre_tabs2, sep="\t", dec=".", quote=F, row.names = F)






cat("WE")

WE4 <- l %>% filter(population=="WE")%>%
  filter(kept==T) %>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>% group_by(memb) %>% summarise(m=mean(value)) %>% arrange(desc(m)) %>% dplyr::select(memb) %>%
  unlist() %>% as.vector() %>% first()

l %>% filter(population=="WE") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",WE4))) %>%
  unlist() %>% as.vector() %>% mean()

l %>% filter(population=="WE") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",WE4))) %>%
  unlist() %>% as.vector() %>% sd()


cat("EE")

EE4 <- l %>% filter(population=="EE")%>%
  filter(kept==T) %>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>% group_by(memb) %>% summarise(m=mean(value)) %>% arrange(desc(m)) %>% dplyr::select(memb) %>%
  unlist() %>% as.vector() %>% first()

l %>% filter(population=="EE") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",EE4))) %>%
  unlist() %>% as.vector() %>% mean()

l %>% filter(population=="EE") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",EE4))) %>%
  unlist() %>% as.vector() %>% sd()


cat("WA")

WA4 <- l %>% filter(population=="WA")%>%
  filter(kept==T) %>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>% group_by(memb) %>% summarise(m=mean(value)) %>% arrange(desc(m)) %>% dplyr::select(memb) %>%
  unlist() %>% as.vector() %>% first()

l %>% filter(population=="WA") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",WA4))) %>%
  unlist() %>% as.vector() %>% mean()

l %>% filter(population=="WA") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",WA4))) %>%
  unlist() %>% as.vector() %>% sd()



cat("EA")

EA4 <- l %>% filter(population=="EA")%>%
  filter(kept==T) %>%
  dplyr::select(ERGE, starts_with("STRmemb4_")) %>%
  dplyr::select((-one_of("STRmemb4_max"))) %>% 
  gather(memb, value, -ERGE) %>% 
  mutate(memb=str_replace_all(memb,"STRmemb4_","")) %>% group_by(memb) %>% summarise(m=mean(value)) %>% arrange(desc(m)) %>% dplyr::select(memb) %>%
  unlist() %>% as.vector() %>% first()

l %>% filter(population=="EA") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",EA4))) %>%
  unlist() %>% as.vector() %>% mean()

l %>% filter(population=="EA") %>%
  filter(kept==T) %>%
  dplyr::select(starts_with(paste0("STRmemb4_",EA4))) %>%
  unlist() %>% as.vector() %>% sd()


sessionInfo()