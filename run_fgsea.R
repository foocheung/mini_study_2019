#!/usr/local/bin/Rscript

library(readr)
library(fgsea)
library(tidyverse)
library(docopt)
library(readxl)
setwd('./')


"
Usage:
  run_fgsea.R [--number=<number>] 
  rand-unif (-h | --help | --version)

Description:   This program generates fgsea
Options:
          --version       Show the current version.
          --number=<num>  [default: 1] Thereshold
          
" -> doc

args    <- docopt(doc)

number <- as.numeric(args $ `--number`)

threshold<-number
dddd<-number
a<-read_delim("./results_ctrl_case.txt",delim=' ')


cat(paste0(runif(number), collapse = '\n'), '\n')


a$fcSign=sign(a$log2FoldChange)
a$logP=-log10(a$pvalue)
a$metric=a$logP/a$fcSign 

#res2<-a %>% arrange(desc(metric)) %>% select(id,metric)resx<-a %>% select(id,stat)
resxx<-a %>% dplyr::select(id,stat)


ranks1<-deframe(as.data.frame(resxx))
tm<- gmtPathways("./BTM_for_GSEA_20131008.gmt")

fgseaRes2 <- fgsea(pathways=tm, stats=ranks1, nperm=10000,maxSize=500)

hu_fgsea<-fgseaRes2 %>% arrange(padj) %>% filter(padj < threshold)

#f<-merge(fgseaRes4,fgseaRes3,by="pathway" )

#ggplot(f,aes(x=f$NES.x,y=f$NES.y)) + geom_point()



ggplot(hu_fgsea, aes(y =reorder(pathway,-log10(padj)) , x = "", size = -log10(padj), fill = NES)) + 
  geom_point(shape = 21) + 
  scale_fill_viridis_c() + xlab("") +
  theme_bw() +
  scale_x_discrete(position = "top") + 
  #theme(axis.text.x=element_text(angle = 0, hjust = 0)) + 
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "right") + 
  labs(size = '-log10(ajdusted P value)', fill = 'Normalized Enrichment Score') +
  theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
  theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
  theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
  guides(shape = guide_legend(override.aes = list(size = 5))) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8))



###############
###############

###############
###############

###############
###############
mouse<-read_tsv("./mouse_data.txt")
mouse<-mouse %>%  mutate(id = toupper(genesymbol))
#resB<-mouse %>% select(id,foldchangeTFAMwt)
mouse$genesymbol<-sub("\\|.*", "\\1", mouse$genesymbol)
resB<-mouse %>% dplyr::select(id,foldchangeTFAMwt)
resB<-as.tibble(resB)

ranksB<-deframe(as.data.frame(resB))
tm<- gmtPathways("./BTM_for_GSEA_20131008.gmt")

fgseaResB <- fgsea(pathways=tm, stats=ranksB, nperm=10000,maxSize=500)


mu_fgsea<-fgseaResB %>% arrange(padj) %>% filter(padj <threshold)




#########################3
mu_fgsea<-mu_fgsea %>% mutate(grp="mouse")
hu_fgsea<-hu_fgsea %>% mutate(grp="human")
allfgsea<-rbind(hu_fgsea,mu_fgsea)


allfgsea$ID<-sub("^.*?\\(([M|S].*)\\)[^)]*$", "\\1", allfgsea$pathway)


##c("Cell Cycle", "0", "T Cells", "NK Cells", 
#"Monocytes", "Inflamatory TLR Chemokines", "Neutrophils", "Interferon/Antiviral Sensing", 
#"DC Activation", "Mitochondria", "B Cells", "Plasma Cells", "Golgi/ER", 
#"Antigen Presentation"))



modules<-read_excel("./mod_comp.xlsx")
BCELLS<-modules %>% filter(modules$`High-level annotation group` == 'B CELLS') %>% dplyr::select(BTM)
allfgsea2 <- allfgsea %>% mutate(Ann = case_when(
  (ID %in% c('M5.0','M95.0','M28','M71','M95.1','M200') ~ "Antigen Presentation"),
  (ID %in% c("M47.0", "M47.1", "M47.2", "M47.3", "M47.4","M54", "M69", "M83", "M9", "M58", "S2", "S8", "S9") ~ "B Cells"),
  (ID %in% c("M156.0", "M156.1", "M123", "S3") ~ "Plasma Cells"),
  (ID %in% c("M61.0", "M61.1", "M61.2", "M7.2", "M157", "S1") ~ "NK Cells"),
  (ID %in% c("M11.0", "M23", "M73", "M81", "M118.0","S4", "M4.15", "M118.1")~ "Monocytes"),
  (ID %in% c(" M4.9", "M7.0", "M7.1", "M7.4", "M14","M19", "M35.0", "M35.1", "M223", "M5.1", "M7.3", "M18", "M36","M44", "M52", "S0", "S6", "S7") ~ "T Cells"),
  (ID %in% c("M4.0", "M4.1", "M4.4", "M4.5", "M4.6","M4.7", "M4.10", "M4.11", "M6", "M8", "M22.0", "M22.1", "M31","M46", "M76", "M103", "M167", "M169", "M230", "M4.8", "M15","M37.3") ~ "Cell Cycle"),
  (ID %in% c("M40", "M43.0", "M43.1", "M53", "M64", "M67", "M86.0", "M86.1", "M119", "M165", "M168") ~  "DC Activation"),
  (ID %in% c("M16", "M25", "M27.0", "M27.1", "M29","M37.0", "M115", "M146", "S5", "S10", "S11", "M33", "M50") ~ "Inflamatory TLR Chemokines"),
  (ID %in% c("M11.2", "M37.1", "M132", "M163") ~ "Neutrophils"),
  (ID %in% c("M75", "M127", "M158.0", "M158.1", "M13", "M68", "M111.0", "M111.1", "M150")~ "Interferon/Antiviral Sensing"),
  (ID %in% c("M187", "M238", "M231","M219","M216")~ "Mitochondria"),
  (ID %in% c("M87","M113","M37.2")~ "Golgi/ER"),
  
  (TRUE  ~ "0")
  
) ) 




allfgsea2$Ann<-factor(allfgsea2$Ann, levels=c("Interferon/Antiviral Sensing","DC Activation","Mitochondria",
                                               "B Cells","Antigen Presentation",
                                               "Cell Cycle", "Golgi/ER", "Monocytes","Neutrophils",
   "T Cells", "NK Cells", "Inflamatory TLR Chemokines",
                        "Plasma Cells" ,"0"
                       ))







p1<-ggplot(allfgsea2, aes(y =reorder(pathway,-log10(padj)) , x = grp, size = -log10(padj), fill = NES)) + 
  geom_point(shape = 21) + 
  scale_fill_viridis_c() + xlab("") +
  theme_bw() +
  scale_x_discrete(position = "top") + 
  #theme(axis.text.x=element_text(angle = 0, hjust = 0)) + 
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "right") + 
  labs(size = '-log10(ajdusted P value)', fill = 'Normalized Enrichment Score') +
  theme(legend.title = element_text(face = "bold",colour = "black", size = 8)) +
  theme(axis.text.y = element_text(size = 8, face = "bold", color = "black")) + 
  theme(axis.text.x = element_text(size = 8.5, face = "bold", color = "black")) + 
  guides(shape = guide_legend(override.aes = list(size = 5))) + 
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) + facet_wrap(Ann~.,scales="free",strip.position="left")

pdf(paste("./mini_tmod_",threshold,".pdf",sep=''),height=18,width=24)
p1
dev.off()




write.table(as.matrix(allfgsea2),paste("./mini_tmod_",threshold,".txt",sep=''),sep="\t")













