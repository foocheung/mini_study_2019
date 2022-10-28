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
a<-read_delim("./results_ctrl_case.txt",delim=' ')


a$fcSign=sign(a$log2FoldChange)
a$logP=-log10(a$pvalue)
a$metric=a$logP/a$fcSign 

#res2<-a %>% arrange(desc(metric)) %>% select(id,metric)resx<-a %>% select(id,stat)
resxx<-a %>% dplyr::select(id,stat)


ranks1<-deframe(as.data.frame(resxx))
tm<- gmtPathways("./c2.cp.kegg.v7.1.symbols.gmt")

fgseaRes2 <- fgsea(pathways=tm, stats=ranks1, nperm=10000,maxSize=500)

hu_fgsea<-fgseaRes2 %>% arrange(padj) %>% filter(padj <threshold)

#f<-merge(fgseaRes4,fgseaRes3,by="pathway" )

#ggplot(f,aes(x=f$NES.x,y=f$NES.y)) + geom_point()





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
tm<- gmtPathways("./c2.cp.kegg.v7.1.symbols.gmt")

fgseaResB <- fgsea(pathways=tm, stats=ranksB, nperm=10000,maxSize=500)


mu_fgsea<-fgseaResB %>% arrange(padj) %>% filter(padj <threshold)




#########################3
mu_fgsea<-mu_fgsea %>% mutate(grp="mouse")
hu_fgsea<-hu_fgsea %>% mutate(grp="human")
allfgsea<-rbind(hu_fgsea,mu_fgsea)


allfgsea$ID<-sub("^.*?\\(([M|S].*)\\)[^)]*$", "\\1", allfgsea$pathway)



##############same_directions

hu_up_fg_id<-allfgsea %>% filter(NES > 0 &  grp == 'human') %>% dplyr::select(ID)
dput(hu_up_fg_id$ID)
allfgsea1<-allfgsea %>% filter(NES > 0 &  grp == 'mouse' & ID %in% dput(hu_up_fg_id$ID))

allfgsea_all1<-allfgsea %>% filter(ID %in% dput(allfgsea1$ID)) %>% mutate(direction="human:up mouse:up")






hu_down_fg_id<-allfgsea %>% filter(NES < 0 &  grp == 'human') %>% dplyr::select(ID)
dput(hu_down_fg_id$ID)
allfgsea2<-allfgsea %>% filter(NES < 0 &  grp == 'mouse' & ID %in% dput(hu_down_fg_id$ID)) 
allfgsea_all2<-allfgsea %>% filter(ID %in% dput(allfgsea2$ID)) %>% mutate(direction="human:down mouse:down")


##############opposites


up_dn_fg_id<-allfgsea %>% filter(NES > 0 &  grp == 'human') %>% dplyr::select(ID)
dput(up_dn_fg_id$ID)
allfgsea3<-allfgsea %>% filter(NES < 0 &  grp == 'mouse' & ID %in% dput(up_dn_fg_id$ID)) 
allfgsea_all3<-allfgsea %>% filter(ID %in% dput(allfgsea3$ID)) %>% mutate(direction="human:up mouse:down")



down_dn_fg_id<-allfgsea %>% filter(NES < 0 &  grp == 'human') %>% dplyr::select(ID)
dput(down_dn_fg_id$ID)
allfgsea4<-allfgsea %>% filter(NES > 0 &  grp == 'mouse' & ID %in% dput(down_dn_fg_id$ID)) 

allfgsea_all4<-allfgsea %>% filter(ID %in% dput(allfgsea4$ID))%>% mutate(direction="human:down mouse:up")


########################################################################

allfgsea_all6<-rbind(allfgsea_all1,allfgsea_all2,allfgsea_all3,allfgsea_all4)



p4<-ggplot(allfgsea_all6, aes(y =reorder(pathway,-log10(padj)) , x = grp, size = -log10(padj), fill = NES)) + 
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
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) + facet_wrap(.~ direction,scales="free",ncol = 2)



pdf(paste("./miniKegg_",threshold,".pdf",sep=''),height=12,width=20)
p4
dev.off()



write.table(as.matrix(allfgsea_all6),paste("./miniKegg_",threshold,".txt",sep=''),sep="\t")






