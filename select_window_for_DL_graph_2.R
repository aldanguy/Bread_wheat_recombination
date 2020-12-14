
cat("\n\nselect_window_for_DL_graph_2.R\n\n")
Sys.time()
rm(list = ls())
set.seed(1)
graphics.off()
variables <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))



titre_ld <- variables[1]
titre_hr <- variables[2]
titre_fenetres <- variables[3]
titre_res <- variables[4]
titre_g <- variables[5]

titre_ld <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/w_WE_3B_86.ld"
titre_hr <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/HR.txt"
titre_fenetres <-  "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/windows.txt"
titre_res <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/PHASE_summary_outputs.txt"
titre_g <- "/home/adanguydesd/Documents/These_Alice/recombinaison/pipeline/020820/graphes/LD.tiff"

hr <- fread(titre_hr)
fr <- fread(titre_fenetres)
ld <- fread(titre_ld)
res <- fread(titre_res)

hr <- hr %>% filter(kept==T) %>%
  filter(population=="WE") %>%
  filter(chr=="3B") %>%
  arrange(desc(lambda_med)) 


w <- fr %>% filter(population=="WE" & chr=="3B") %>%
  inner_join(hr, by="chr") %>%
  filter(posl >= centerlfen & posr <= centerrfen) %>%
  arrange(desc(nb_mqs), desc(lambda_med)) %>%
  slice(1) %>%
  dplyr::select(ID) %>%
  unlist %>%
  as.vector()


region <- fr %>% filter(population=="WE" & chr=="3B") %>%
  inner_join(hr, by="chr") %>%
  filter(posl >= centerlfen & posr <= centerrfen) %>%
  arrange(desc(nb_mqs), desc(lambda_med)) %>%
  slice(1) %>%
  dplyr::select(region.x)%>%
  unlist %>%
  as.vector()

res <- res %>% filter(ID==w)

poshr <- fr %>% filter(ID %in% w) %>% inner_join(hr) %>%
  filter(posl >= centerlfen & posr <= centerrfen) %>%
  mutate(chr=8) %>%
  dplyr::select(chr, posl, posr, lambda_med)


order <- ld %>% pivot_longer(cols = c("SNP_A","SNP_B")) %>% filter(R2==1) %>%
  mutate(pos=ifelse(name=="SNP_A",BP_A,BP_B)) %>%
  dplyr::select(value, pos) %>% unique() %>%
  arrange(pos) %>% 
  mutate(ordre=c(1:nrow(.))) 



ld2 <- ld %>% arrange(BP_A, BP_B) %>%
  full_join(poshr, by=c("CHR_A"="chr")) %>%
  mutate(hr=ifelse(BP_A == posl & BP_B==posr,T,NA)) %>%
  dplyr::select(SNP_A,BP_A, BP_B,SNP_B,R2, hr,lambda_med) %>%
  arrange(BP_A, BP_B) %>%
  mutate(SNP_A=factor(SNP_A, levels=order$value)) %>%
  mutate(SNP_B=factor(SNP_B, levels=order$value)) %>%
  mutate(lambda_med=round(lambda_med)) %>%
  mutate(lambda_med=factor(lambda_med, levels=rev(sort(unique(round(poshr$lambda_med))))))
  

ga <- ld2 %>%
  ggplot(aes(x=SNP_A, y=SNP_B, fill=R2)) +
  theme_void() +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "blue", limit = c(0,1), space = "Lab",
                       name="R2")+
  geom_point(data=ld2 %>% filter(hr==T), aes(x=SNP_A, y=SNP_B, col=lambda_med), size=5)+
  geom_point(data=ld2 %>% filter(hr==T), aes(x=SNP_B, y=SNP_A,  col=lambda_med), size=5)+
  guides(colour=guide_legend(expression(lambda))) +
  ggtitle(paste0("SNP window ", w, " (",region,")"))+
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16), 
        legend.position = "left",plot.margin = unit(c(0,1,1,0), "cm"))


gb <- res %>% ggplot(aes(x=posSNPlpop, y=lambda_rho_med)) + geom_line() +
  theme_light() +
  scale_x_continuous(labels=function(x)x/1e6) +
  scale_y_continuous(labels=function(x)x*1e3) +
  ylab(expression(paste(lambda, " * ",rho, " (/kb)"))) +
  xlab("physical position (Mb)") +
  geom_point(data=ld2 %>% filter(hr==T), aes(x=BP_A, y=0,col=lambda_med), size=5)+
  guides(colour=guide_legend(expression(lambda))) +
  ggtitle(paste0("LD-based recombination landscape "))+
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.title.x = element_text(size=18))


g <- ggarrange(ga,gb, ncol=2, nrow=1, labels=c("A","B"), font.label=list(size=18))
g
tiff(titre_g, units="in", width = 12, height = 6, res=200, compression = "lzw")
g
dev.off()



sessionInfo()
