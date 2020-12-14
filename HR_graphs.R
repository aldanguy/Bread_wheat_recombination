







Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(viridis))


cat("\n HR_graphs.R\n")
cat("\n\nVariables : \n")
print(variables)
cat("\n\n")


titre_coloc_permutation <- variables[1]
titre_HR <- variables[2]
titre_lambda <-variables[3]
titre_coloc_ref <-variables[4]
titre_resume <- variables[5]
titre_clique_HR_1000 <-variables[6]
titre_graphe_coloc <- variables[7]
titre_graphe_distances <- variables[8]
titre_f5 <- variables[9]
titre_graphe_relativ2 <- variables[10]

titre_sortie_temp1 <-variables[11]
titre_sortie_temp2 <-variables[12]

titre_graph_relativ3 <- variables[13]



    


# titre_coloc_permutation <- "/work/adanguy/these/pipeline/020820/HR/permutations.txt"                             
# titre_HR <- "/work/adanguy/these/pipeline/020820/HR/cliques_HR.txt"                               
# titre_lambda <- "/work/adanguy/these/pipeline/020820/HR/HR.txt"                                       
# titre_coloc_ref <- "/work/adanguy/these/pipeline/020820/HR/coloc_HR.txt"                                 
# titre_resume <- "/work/adanguy/these/pipeline/020820/PHASE/PHASE_summary_outputs.txt"                 
# titre_clique_HR_1000 <- "/work/adanguy/these/pipeline/020820/HR/cliques_HR_1000.txt"                          
# titre_graphe_coloc <- "/work/adanguy/these/pipeline/020820/graphs/coloc_HR.png"                             
# titre_graphe_distances <- "/work/adanguy/these/pipeline/020820/graphs/coloc_HR_distances.png"                   
# titre_f5 <-"/work/adanguy/these/pipeline/020820/graphs/relative_lambda_HR_populations.png"       
# titre_graphe_relativ2 <-  "/work/adanguy/these/pipeline/020820/graphs/relative_lambda_HR_populations_2.png"     
# titre_sortie_temp1 <- "/work/adanguy/these/pipeline/020820/HR/relativ_lambda_HR_populations.txt"            
# titre_sortie_temp2 <-  "/work/adanguy/these/pipeline/020820/HR/relativ_lambda_HR_combination_populations.txt"
# titre_graph_relativ3 <-"/work/adanguy/these/pipeline/020820/graphs/relative_lambda_HR_populations_3.png"  

  
  #  titre_sortie_temp1 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/relativ_lambda_HR_populations.txt"
  #  titre_f5 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/Figure6.tiff"
  #  titre_graphe_relativ2  <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/relative_lambda_HR_2.png"
  #   
  #   
  #   
  #  titre_coloc_permutation <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/permutations.txt"
  #  titre_HR <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/cliques_HR.txt"
  # titre_coloc_ref <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/coloc_HR.txt" 
  #  titre_lambda <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/HR.txt"
  #  titre_resume <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/PHASE_summary_outputs.txt"
  #  titre_clique_HR_1000 <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/cliques_HR_1000.txt"    

 
 distance_max_HR <- 1e5
 step <- 1e3
 
cat("\n\n Input 1 : permutation of HRs \n\n")
permutation <- fread(titre_coloc_permutation) 
head(permutation)
permutation <- permutation %>%
  mutate(v=n/tcliques)


permutation %>% group_by(coloc, iteration) %>%
  filter(coloc %in% c("WE","EE","WA","EA")) %>%
  mutate(prop_coloc=n/tcliques) %>%
  group_by(iteration) %>%
  summarise(prop_coloc=sum(prop_coloc)) %>%
  ungroup() %>%
  summarise(mean_prop_specific=mean(prop_coloc), sd_prop_specific=sd(prop_coloc))


permutation %>% group_by(coloc, iteration) %>%
  filter(coloc %in% c("WEEEWAEA")) %>%
  mutate(prop_coloc=n/tcliques) %>%
  group_by(iteration) %>%
  summarise(prop_coloc=sum(prop_coloc)) %>%
  ungroup() %>%
  summarise(mean_prop_commun4=mean(prop_coloc), sd_prop_commun4=sd(prop_coloc))


cat("\n\n Input 2 : coloc of true HRs \n\n")
coloc <- fread(titre_coloc_ref) 
coloc
coloc <- coloc  %>%
  mutate(vref=n/tcliques)





cat("\n\n Input 3 : coloc of HR_1000 \n\n")
sh <- fread(titre_clique_HR_1000) 
head(sh)



# cat("\n\n Input 4 : values of lambda for HR \n\n")
lambda <- fread(titre_lambda) 
head(lambda)


cat("\n\n Input 5 : HR positions \n\n")
HR <- fread(titre_HR) 
head(HR)
# lambda <- HR



cat("\n\n Input 6 : summary outputs from PHASE \n\n")
res <- fread(titre_resume) 
head(res)




color=viridis(5)


coloc <- coloc %>% mutate(color=case_when(coloc %in% c("WEEEWAEA")~color[4],
                                          coloc %in% c("WEEEWA","WEEEEA","EEWAEA","WEWAEA")~color[3],
                                          coloc %in% c("WEEE","WEWA","WEEA","EEWA","EEEA","WAEA")~color[2],
                                          coloc %in% c("WE","EE","WA","EA")~color[1]
)) %>%
  mutate(coloc=factor(coloc, levels=c("WEEEWAEA",
                                      "WEEEWA","WEEEEA","WEWAEA","EEWAEA",
                                      "WEEE","WEWA","WEEA","EEWA","EEEA","WAEA",
                                      "WE","EE","WA","EA"))) %>%
  arrange(coloc)
# 
# test <- coloc$color
# names(test) <- coloc$coloc

p1 <- upset(sh, keep.order=T,
            sets = rev(c("WE","EE","WA","EA")), order.by = c("degree" ),
            sets.bar.color=c("#2B83BA","#ABDDA4","#FDAE61","#D7191C"),
            mainbar.y.label = "Number of shared HR",
            queries = list(list(query = elements, params = list("WE"), color =color[1], active = TRUE),
                           list(query = elements, params = list("EE"), color =color[1], active = TRUE),
                           list(query = elements, params = list("WA"), color =color[1], active = TRUE),
                           list(query = elements, params = list("EA"), color =color[1], active = TRUE),
                           list(query = elements, params = list("WE","EE"), color =color[2], active = TRUE),
                           list(query = elements, params = list("WE","WA"), color =color[2], active = TRUE),
                           list(query = elements, params = list("WE","EA"), color =color[2], active = TRUE),
                           list(query = elements, params = list("EE","WA"), color =color[2], active = TRUE),
                           list(query = elements, params = list("EE","EA"), color =color[2], active = TRUE),
                           list(query = elements, params = list("WA","EA"), color =color[2], active = TRUE),
                           list(query = elements, params = list("WE","EE","WA"), color =color[3], active = TRUE),
                           list(query = elements, params = list("WE","EE","EA"), color =color[3], active = TRUE),
                           list(query = elements, params = list("WE","WA","EA"), color =color[3], active = TRUE),
                           list(query = elements, params = list("EE","WA","EA"), color =color[3], active = TRUE),
                           list(query = elements, params = list("WE","EE","WA","EA"), color =color[4], active = TRUE)
            ))



p1 <- cowplot::plot_grid(NULL, p1$Main_bar, p1$Sizes, p1$Matrix,
                         nrow=2, align='hv', rel_heights = c(3,1),
                         rel_widths = c(2,3))


table(sh$coloc)
# p1
# png(titre_graphe)
# p1
# dev.off()

#ABDDA4 = green
#D7191C = blue
#FDAE61 = red
#2B83BA = orange



p2 <- permutation %>%
  full_join(coloc, by=c("coloc")) %>%
  mutate(coloc=factor(coloc, levels=c("WEEEWAEA",
                                      "WEEEWA","WEEEEA","WEWAEA","EEWAEA",
                                      "WEEE","WEWA","EEWA","WEEA","EEEA","WAEA",
                                      "WE","EE","WA","EA"))) %>%
  ggplot(aes(x=coloc, y =v)) + geom_boxplot(col="grey")  +
  geom_point(aes(x=coloc, y=vref, col=coloc), shape=18, size=5, show.legend = F) +
  scale_color_manual(values=c("WEEEWAEA"=color[4],
                              "WEEEWA"=color[3],
                              "WEEEEA"=color[3],
                              "WEWAEA"=color[3],
                              "EEWAEA"=color[3],
                              "WEEE"=color[2],
                              "WEWA"=color[2],
                              "EEWA"=color[2],
                              "WEEA"=color[2],
                              "EEEA"=color[2],
                              "WAEA"=color[2],
                              "WE"=color[1],
                              "EE"=color[1],
                              "WA"=color[1],
                              "EA"=color[1]                              )) +
  ylab("Expected by chance (%)") +
  xlab("")  +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 




# f5a <- permutation %>%
#   full_join(coloc, by=c("coloc")) %>%
#   mutate(coloc=case_when(coloc=="WEEE"~"WE-EE",
#                    coloc=="WEWA"~"WE-WA",
#                    coloc=="WEEA"~"WE-EA",
#                    coloc=="EEWA"~"EE-WA",
#                    coloc=="EEEA"~"EE-EA",
#                    coloc=="WAEA"~"WA-EA",
#                    coloc=="WEEEWA"~"WE-EE-WA",
#                    coloc=="WEEEEA"~"WE-EE-EA",
#                    coloc=="WEWAEA"~"WE-WA-EA",
#                    coloc=="EEWAEA"~"EE-WA-EA",
#                    
#                    coloc=="WE"~"WE",
#                    coloc=="EE"~"EE",
#                    coloc=="WA"~"WA",
#                    coloc=="EA"~"EA",
#                    
#                    coloc=="WEEEWAEA"~"WE-EE-WA-EA"))%>%
#   mutate(coloc=factor(coloc, levels=c(c("WE-EE-WA-EA",
#                                       rev(c("EE-WA-EA","WE-EE-EA","WE-WA-EA","WE-EE-WA")),
#                                       rev(c("EE-EA","WE-EA", "WA-EA","EE-WA","WE-EE","WE-WA")),
#                                      rev(c("EA","EE","WE","WA")))))) %>%
#   mutate(groupe=case_when(coloc%in%c("WE","EE","WA","EA")~"population specific",
#                           coloc %in% c( "WE-EE","WE-WA","WE-EA","EE-WA","EE-EA","WA-EA") ~ "2 populations" ,
#                           coloc %in% c("WE-EE-WA","WE-EE-EA","WE-WA-EA","EE-WA-EA")~"3 populations",
#                           coloc %in% c("WE-EE-WA-EA")~ "4 populations"))%>%
#   mutate(groupe=factor(groupe, levels = rev(c("4 populations","3 populations","2 populations","population specific")))) %>%
#   mutate(temp="temp") %>%
#   ggplot(aes(x=coloc, y =v, fill=temp)) + geom_boxplot(col="grey", show.legend = T)  +
#   scale_fill_manual(values="grey", name="", labels="by chance")+
#   
#   geom_point(aes(x=coloc, y=vref, col=coloc), shape=18, size=5, show.legend = F) +
#   scale_color_manual(values=c("WE-EE-WA-EA"=color[4],
#                               "WE-EE-WA"=color[3],
#                               "WE-EE-EA"=color[3],
#                               "WE-WA-EA"=color[3],
#                               "EE-WA-EA"=color[3],
#                               "WE-EE"=color[2],
#                               "WE-WA"=color[2],
#                               "EE-WA"=color[2],
#                               "WE-EA"=color[2],
#                               "EE-EA"=color[2],
#                               "WA-EA"=color[2],
#                               "WE"=color[1],
#                               "EE"=color[1],
#                               "WA"=color[1],
#                               "EA"=color[1]                              ), guide=F) +
#   ylab("Percentage of HR shared by") +
#     facet_grid(.~groupe, scales="free")+
#   xlab("")  +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
#   geom_text(aes(x=coloc, y=vref, col=coloc, label=round(vref, digits=2)),hjust=0.8, vjust=-1, check_overlap = TRUE, size=5) +
#   theme(legend.position = "right",
#           strip.text.y = element_text(size = 15, color="black"),
#                 strip.text.x = element_text(size = 20, color="black"),
#                 strip.background = element_rect(color="black", fill="white"),
#                 axis.title.y = element_text(size=22),
#         axis.text.x = element_text(size=20),
#         axis.text.y = element_text(size=20),
#                 axis.title.x = element_text(size=22),
#                 legend.text=element_text(size=20),
#                 legend.title=element_text(size=22),
#         panel.grid.minor = element_blank())



data <- permutation %>%
  full_join(coloc, by=c("coloc")) %>%
  mutate(coloc=case_when(coloc=="WEEE"~"WE-EE",
                         coloc=="WEWA"~"WE-WA",
                         coloc=="WEEA"~"WE-EA",
                         coloc=="EEWA"~"EE-WA",
                         coloc=="EEEA"~"EE-EA",
                         coloc=="WAEA"~"WA-EA",
                         coloc=="WEEEWA"~"WE-EE-WA",
                         coloc=="WEEEEA"~"WE-EE-EA",
                         coloc=="WEWAEA"~"WE-WA-EA",
                         coloc=="EEWAEA"~"EE-WA-EA",
                         
                         coloc=="WE"~"WE",
                         coloc=="EE"~"EE",
                         coloc=="WA"~"WA",
                         coloc=="EA"~"EA",
                         
                         coloc=="WEEEWAEA"~"WE-EE-WA-EA"))%>%
  mutate(coloc=factor(coloc, levels=c(c("WE-EE-WA-EA",
                                        rev(c("EE-WA-EA","WE-EE-EA","WE-WA-EA","WE-EE-WA")),
                                        rev(c("EE-EA","WE-EA", "WA-EA","EE-WA","WE-EE","WE-WA")),
                                        rev(c("EA","EE","WE","WA")))))) %>%
  mutate(groupe=case_when(coloc%in%c("WE","EE","WA","EA")~"population specific",
                          coloc %in% c( "WE-EE","WE-WA","WE-EA","EE-WA","EE-EA","WA-EA") ~ "2 populations" ,
                          coloc %in% c("WE-EE-WA","WE-EE-EA","WE-WA-EA","EE-WA-EA")~"3 populations",
                          coloc %in% c("WE-EE-WA-EA")~ "4 populations"))%>%
  mutate(groupe=factor(groupe, levels = rev(c("4 populations","3 populations","2 populations","population specific")))) %>%
  group_by(groupe, iteration) %>%
  mutate(vref=sum(vref), v=sum(v)) %>%
  mutate(temp="temp") 

data_pop_spe <- data %>% filter(groupe%in%"population specific") %>% dplyr::select(groupe, iteration, v, vref) %>%
  as.data.frame()
data_pop_2 <- data %>% filter(!groupe%in%"population specific")  %>%
  dplyr::select(groupe, iteration, v, vref) %>% unique() %>% 
  arrange(iteration, groupe) %>%
  group_by(iteration) %>%
  summarise(v=sum(v), vref=sum(vref)) %>%
  mutate(groupe="2 populations")%>% 
  dplyr::select(groupe, iteration, v, vref)%>%
  as.data.frame()
data_pop_3 <- data %>% filter(!groupe %in%c("population specific","2 populations"))  %>%
  dplyr::select(groupe, iteration, v, vref) %>% unique() %>% 
  arrange(iteration, groupe) %>%
  group_by(iteration) %>%
  summarise(v=sum(v), vref=sum(vref)) %>%
  mutate(groupe="3 populations")%>% 
  dplyr::select(groupe, iteration, v, vref) %>%
  ungroup()%>%
  as.data.frame()
data_pop_4 <- data %>% filter(groupe%in%"4 populations") %>% dplyr::select(groupe, iteration, v, vref)%>%
  as.data.frame()
data <- rbind(data_pop_spe, data_pop_2, data_pop_3, data_pop_4) %>%
  mutate(temp="temp") %>%
  mutate(groupe=factor(groupe, levels = rev(c("4 populations","3 populations","2 populations","population specific"))))



f5a <- data %>% ggplot(aes(x=groupe, y =v, fill=temp)) + geom_boxplot(col="grey", show.legend = T)  +
  scale_fill_manual(values="grey", name="", labels="by chance")+
  
  geom_point(aes(x=groupe, y=vref, col=groupe), shape=18, size=5, show.legend = F) +
  scale_color_manual(values=c("4 populations"=color[4],
                              "3 populations"=color[3],
                              "2 populations"=color[2],
                              "population specific"=color[1]                              ), guide=F) +
  ylab("Percentage of HRIs shared by") +
  facet_grid(.~groupe, scales="free")+
  xlab("")  +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  geom_text(aes(x=groupe, y=vref, col=groupe, label=round(vref, digits=2)),hjust=0.8, vjust=-1, check_overlap = TRUE, size=10) +
  theme(legend.position = "right",
        strip.text.y = element_text(size = 15, color="black"),
        strip.text.x = element_text(size = 26, color="black"),
        strip.background = element_rect(color="black", fill="white"),
        axis.title.y = element_text(size=26),
        axis.text.x = element_text(size=0),
        axis.text.y = element_text(size=24),
        axis.title.x = element_text(size=22),
        legend.text=element_text(size=26),
        legend.title=element_text(size=22),
        panel.grid.minor = element_blank())
f5a

cat("\n\n quantile \n\n")
permutation %>%
  full_join(coloc, by=c("coloc")) %>%
  mutate(coloc=factor(coloc, levels=c("WEEEWAEA",
                                      "WEEEWA","WEEEEA","WEWAEA","EEWAEA",
                                      "WEEE","WEWA","WEEA","EEWA","EEEA","WAEA",
                                      "WE","EE","WA","EA"))) %>%
  group_by(coloc) %>%
  summarise(q=ecdf(v)(unique(vref))) %>%
  as.data.frame()



##################################### corriger pour répétitions


p3 <-  HR %>%
  mutate(coloc=factor(coloc, levels=c("WEEEWAEA",
                                      "WEEEWA","WEEEEA","WEWAEA","EEWAEA",
                                      "WEEE","WEWA","EEWA","WEEA","EEEA","WAEA",
                                      "WE","EE","WA","EA"))) %>%
  ggplot(aes(y=lambda_med, x=coloc, col=coloc)) + geom_boxplot(show.legend = F) +
  scale_color_manual(values=c("WEEEWAEA"=color[4],
                              "WEEEWA"=color[3],
                              "WEEEEA"=color[3],
                              "WEWAEA"=color[3],
                              "EEWAEA"=color[3],
                              "WEEE"=color[2],
                              "WEWA"=color[2],
                              "EEWA"=color[2],
                              "WEEA"=color[2],
                              "EEEA"=color[2],
                              "WAEA"=color[2],
                              "WE"=color[1],
                              "EE"=color[1],
                              "WA"=color[1],
                              "EA"=color[1]                              )) +
  scale_y_continuous(trans="log10") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") +
  ylab(expression(lambda)) 



lambda <- HR %>% mutate(nbpop=case_when(coloc=="WEEEWAEA"~4,
                                            coloc %in% c("WEEEWA","WEEEEA","WEWAEA","EEWAEA")~3,
                                            coloc %in% c("WEEE","WEWA","WEEA","EEWA","EEEA","WAEA")~2,
                                            coloc %in% c("WE","EE","WA","EA")~1))

pairwise.wilcox.test(lambda$lambda_med, lambda$nbpop, p.adjust.method = "bonf")


lambda %>% group_by(nbpop) %>% summarise(mediane=median(lambda_med))




p4 <-  lambda %>%
  mutate(coloc=factor(coloc, levels=c("WEEEWAEA",
                                      "WEEEWA","WEEEEA","WEWAEA","EEWAEA",
                                      "WEEE","WEWA","EEWA","WEEA","EEEA","WAEA",
                                      "WE","EE","WA","EA"))) %>%
  ggplot(aes(y=d_genet_histo, x=coloc, col=coloc)) + geom_boxplot(show.legend = F) +
  scale_color_manual(values=c("WEEEWAEA"=color[4],
                              "WEEEWA"=color[3],
                              "WEEEEA"=color[3],
                              "WEWAEA"=color[3],
                              "EEWAEA"=color[3],
                              "WEEE"=color[2],
                              "WEWA"=color[2],
                              "EEWA"=color[2],
                              "WEEA"=color[2],
                              "EEEA"=color[2],
                              "WAEA"=color[2],
                              "WE"=color[1],
                              "EE"=color[1],
                              "WA"=color[1],
                              "EA"=color[1]                              )) +
  scale_y_continuous(trans="log10") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") +
  ylab(expression(paste(lambda,"*", rho,"*d"))) 


p5 <-  lambda %>%
  mutate(coloc=factor(coloc, levels=c("WEEEWAEA",
                                      "WEEEWA","WEEEEA","WEWAEA","EEWAEA",
                                      "WEEE","WEWA","EEWA","WEEA","EEEA","WAEA",
                                      "WE","EE","WA","EA"))) %>%
  ggplot(aes(y=taille, x=coloc, col=coloc)) +geom_boxplot(show.legend = F) +
  scale_color_manual(values=c("WEEEWAEA"=color[4],
                              "WEEEWA"=color[3],
                              "WEEEEA"=color[3],
                              "WEWAEA"=color[3],
                              "EEWAEA"=color[3],
                              "WEEE"=color[2],
                              "WEWA"=color[2],
                              "EEWA"=color[2],
                              "WEEA"=color[2],
                              "EEEA"=color[2],
                              "WAEA"=color[2],
                              "WE"=color[1],
                              "EE"=color[1],
                              "WA"=color[1],
                              "EA"=color[1]                              )) +
  scale_y_continuous(trans="log10") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") +
  ylab("d (pb)") 

#graphe <- ggarrange(uu_c, p2, p3, ncol=2, labels = c("A","B"))


graphe_coloc <- ggarrange(p1,                                                 # First row with scatter plot
                          ggarrange(p2, p3, ncol = 2, labels = c("C", "D")), # Second row with box and dot plots
                          nrow = 2, 
                          labels = c("A","B")                                      # Labels of the scatter plot
) 

graphe_coloc

graphe_distances <- ggarrange(p4,p5, nrow = 2, labels=c("A","B"))
graphe_distances

cat("\n\n Graph of HR coloc \n\n ")
ggsave(titre_graphe_coloc,graphe_coloc, width = 20, height = 20, units = c("cm"))


cat("\n\n Graph of HR coloc particularities \n\n ")
ggsave(titre_graphe_distances,graphe_distances, width = 20, height = 20, units = c("cm"))







res <- res %>%
  na.omit() %>%
  filter(w_center==T) %>%
  dplyr::select(population, chr, posSNPlpop, posSNPrpop, lambda_med) %>%
  mutate(posSNPlpop=posSNPlpop+1) %>%
  pivot_longer(cols=c(posSNPlpop, posSNPrpop), values_to = "pos") %>%
  arrange(population, chr, pos) %>%
  dplyr::select(population, chr, pos, lambda_med) %>%
  filter(! row_number() %in% which(duplicated(.)))

HR <- HR %>% filter(kept==T) %>%
  mutate(pos=round((posr+posl)/2,0)) %>%
  dplyr::select(population, chr, pos, lambda_med) %>%
  arrange(population, chr,pos)


head(HR)
head(res)

chr="3A"
k=0
for (chr in unique(HR$chr)) {
  
  
  resWE <- res %>% filter(population=="WE" & chr==!!chr)
  resEE <- res %>% filter(population=="EE" & chr==!!chr)
  resWA <- res %>% filter(population=="WA" & chr==!!chr)
  resEA <- res %>% filter(population=="EA" & chr==!!chr)
  
  interpolationWE <- approxfun(x=resWE$pos,y=resWE$lambda_med)
  interpolationEE <- approxfun(x=resEE$pos,y=resEE$lambda_med)
  interpolationWA <- approxfun(x=resWA$pos,y=resWA$lambda_med)
  interpolationEA <- approxfun(x=resEA$pos,y=resEA$lambda_med)
  
  HR2 <- HR %>% filter(chr==!!chr)
  
  
  
  
  
  
  
  # plot(phy,a, type = "l")
  
  i=1
  for (i in 1:nrow(HR2)) {
    
    k=k+1
    
    
    population_reference <- HR2[i,"population"]
    
    pos_reference <- HR2[i,"pos"]
    
    lambda_reference <- HR2[i,"lambda_med"]
    
    minmax_reference <- HR2 %>%
      slice(i) %>%
      summarise(posmin=min(pos) - distance_max_HR,posmax=max(pos)+distance_max_HR)
    
    phy <- seq(minmax_reference$posmin,minmax_reference$posmax,step)
    
    n <- length(phy)
    
    aWE <- interpolationWE(phy)/lambda_reference
    aEE <- interpolationEE(phy)/lambda_reference
    aWA <- interpolationWA(phy)/lambda_reference
    aEA <- interpolationEA(phy)/lambda_reference
    
    
    tab <- data.frame(population_ref=rep(population_reference, times=4*n),
                      population=rep(c("WE","EE","WA",'EA'), each=n),
                      chr=rep(chr, times=4*n),
                      pos= rep(phy - pos_reference , times=4),
                      rel_lambda= c(aWE,aEE,aWA,aEA)) %>%
      na.omit()
    

    
    if( k==1) {
      
      cat("\n\n Interpolation of lambda around HR \n")
      print(head(tab))
      write.table(tab, titre_sortie_temp1, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
      
    } else{
      
      write.table(tab, titre_sortie_temp1, col.names = F, row.names = F, dec=".", sep="\t", quote=F, append = T)
      
      
    }
    
  }
  
}


tab <- fread(titre_sortie_temp1)

cat("\n\n Increase of intensity nearby hotspot center \n\n")
tab %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
  mutate(population_ref=factor(population_ref, levels=c("WE","EE","WA","EA"))) %>%
  filter(population_ref != population) %>%
  filter(pos==0 | pos == 1e5 | pos==-1e5) %>%
  mutate(pos=ifelse(pos==-1e5, 1e5, pos))  %>%
  group_by(pos) %>%
  summarise(rel_lambda=median(rel_lambda)) %>%
  mutate(increase=round(rel_lambda-lead(rel_lambda),2))

tab %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
  mutate(population_ref=factor(population_ref, levels=c("WE","EE","WA","EA"))) %>%
  filter(population_ref != population) %>%
  filter(pos==0 | pos == 1e5 | pos==-1e5) %>%
  mutate(pos=ifelse(pos==-1e5, 1e5, pos))  %>%
  group_by(pos, population_ref, population) %>%
  filter(population_ref=="WE") %>%
  summarise(rel_lambda=median(rel_lambda))

f5b <- tab %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
  mutate(population_ref=factor(population_ref, levels=c("WE","EE","WA","EA"))) %>%
  group_by(population, population_ref, pos) %>%
  summarise(rel_lambda=median(rel_lambda)) %>%
  ggplot(aes(x=pos, y=rel_lambda, col=population)) +
  geom_line(position=position_dodge(width=1000), size=2) +
  facet_wrap(population_ref~., ncol=4, scales="free") +
  theme_light() +
  xlab("distance to HRI center (kb)") +
  ylab(expression(paste(lambda, " tested pop / ",lambda," reference pop (%)"))) +
  scale_colour_manual(values=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA"), name="tested\npopulation") +
  scale_x_continuous(breaks = seq(-1e5, 1e5, 1e5), labels=function(x)x/1e3)  +
  scale_y_continuous(labels=function(x)x*100, limits=c(0.11,1), trans="sqrt", breaks = c(0,0.25,0.5,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0,0,1,0), "cm"),
        strip.text.y = element_text(size = 20, color="black"),
                strip.text.x = element_text(size = 26, color="black"),
                strip.background = element_rect(color="black", fill="white"),
                axis.title.y = element_text(size=26),
                axis.title.x = element_text(size=26),
        axis.text.y = element_text(size=24),
        axis.text.x = element_text(size=24),
                legend.text=element_text(size=26),
                legend.title=element_text(size=26),
        panel.grid.minor = element_blank()) 


f5 <- ggarrange(f5a,f5b, nrow=2, labels = c("A","B"), font.label = list(size=30))
f5


tiff(titre_f5, units="in", width = 18, height = 18, res=200, compression = "lzw")
f5
dev.off()



#ABDDA4 = green
#D7191C = blue
#FDAE61 = red
#2B83BA = orange


HR <- lambda %>% filter(kept==T) %>%
  mutate(pos=round((posr+posl)/2,0)) %>%
  group_by(clique) %>%
  mutate(lambda_med=median(lambda_med)) %>%
  ungroup() %>%
  dplyr::select(coloc, chr, pos, lambda_med) %>%
  unique() %>%
  arrange(coloc, chr,pos) %>%
  rename(population=coloc) %>%
  as.data.frame()


head(HR)
head(res)

chr="3A"
k=0
for (chr in unique(HR$chr)) {
  
  
  resWE <- res %>% filter(population=="WE" & chr==!!chr)
  resEE <- res %>% filter(population=="EE" & chr==!!chr)
  resWA <- res %>% filter(population=="WA" & chr==!!chr)
  resEA <- res %>% filter(population=="EA" & chr==!!chr)
  
  interpolationWE <- approxfun(x=resWE$pos,y=resWE$lambda_med)
  interpolationEE <- approxfun(x=resEE$pos,y=resEE$lambda_med)
  interpolationWA <- approxfun(x=resWA$pos,y=resWA$lambda_med)
  interpolationEA <- approxfun(x=resEA$pos,y=resEA$lambda_med)
  
  HR2 <- HR %>% filter(chr==!!chr)
  
  
  
  
  
  
  
  # plot(phy,a, type = "l")
  
  i=1
  for (i in 1:nrow(HR2)) {
    
    k=k+1
    
    
    population_reference <- HR2[i,"population"]
    
    pos_reference <- HR2[i,"pos"]
    
    lambda_reference <- HR2[i,"lambda_med"]
    
    minmax_reference <- HR2 %>%
      slice(i) %>%
      summarise(posmin=min(pos) - distance_max_HR,posmax=max(pos)+distance_max_HR)
    
    phy <- seq(minmax_reference$posmin,minmax_reference$posmax,step)
    
    n <- length(phy)
    
    aWE <- interpolationWE(phy)/lambda_reference
    aEE <- interpolationEE(phy)/lambda_reference
    aWA <- interpolationWA(phy)/lambda_reference
    aEA <- interpolationEA(phy)/lambda_reference
    
    
    tab <- data.frame(population_ref=rep(population_reference, times=4*n),
                      population=rep(c("WE","EE","WA",'EA'), each=n),
                      chr=rep(chr, times=4*n),
                      pos= rep(phy - pos_reference , times=4),
                      rel_lambda= c(aWE,aEE,aWA,aEA)) %>%
      na.omit()

    
    if( k==1) {
      
      cat("\n\n Interpolation of lambda around HR \n")
      print(head(tab))
      write.table(tab, titre_sortie_temp2, col.names = T, row.names = F, dec=".", sep="\t", quote=F)
      
    } else{
      
      write.table(tab, titre_sortie_temp2, col.names = F, row.names = F, dec=".", sep="\t", quote=F, append = T)
      
      
    }
    
  }
  
}


tab <- fread(titre_sortie_temp2)


p7 <- tab %>%
  mutate(population=factor(population, levels=c("WE","EE","WA","EA"))) %>%
  mutate(population_ref=factor(population_ref, levels=c("WE","EE","WA","EA",
                                                        "WEEE","WEWA","WEEA","EEWA","EEEA", "WAEA",
                                                        "WEEEWA","WEEEEA","WEWAEA","EEWAEA",
                                                        "WEEEWAEA"))) %>%
  group_by(population, population_ref, pos) %>%
  summarise(rel_lambda=mean(rel_lambda)) %>%
  ggplot(aes(x=pos, y=rel_lambda, col=population)) +
  geom_line() +
  facet_wrap(population_ref~., ncol=15) +
  theme_light() +
  xlab("distance to clique center (pb)") +
  ylab(expression(paste("relative ", lambda))) +
  scale_colour_manual(values=c("#D7191C","#FDAE61","#ABDDA4","#2B83BA")) +
  theme(legend.position = "bottom")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_continuous(breaks = seq(-1e5, 1e5, 1e5))



cat("\n\n Graph relativ lambda to clique per population \n")
ggsave(titre_graphe_relativ2, p7, width = 15)



sessionInfo()
