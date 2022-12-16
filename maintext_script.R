library(tidyverse)
library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)
library(biomaRt)
library(EnhancedVolcano)
library(tibble)
library(limma)
library(tximport)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tmod)
library(fgsea)
library(org.Hs.eg.db)
library(viridis)
library(clusterProfiler)
library(GOplot)

setwd("~/Desktop/MINIfinal")

###Prep for importing files from RSEM counts 

samples <- read.table(file.path("test names.txt"), header = T)
files <- file.path("~/Desktop/RSEM_results_copy", samples$file.name)
names(files) <-paste0("sample", 1:81)
all(file.exists(files)) ##TRUE

txi.rsem <- tximport(files, type = "rsem", txIn = F, txOut = F)
txi.rsem$length[txi.rsem$length == 0] <- 1

patient.meta<-read.csv("CHI sample IDs.csv")
sampleTable <- patient.meta[,3:9]
rownames(sampleTable) <- colnames(txi.rsem$counts)

all(rownames(sampleTable) == colnames(txi.rsem$counts)) ##TRUE

##Create DESeq object
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~Group)
dds <- estimateSizeFactors(dds)
dds$Group<-factor(dds$Group, levels = c("Control", "Case"))
dds<-dds[rowSums(counts(dds))>=2,]

vsd<-vst(dds, blind=T)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
df<- as.data.frame(colData(dds)[,c("Group", "Gender")])
pheatmap<-pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, 
         annotation_col=df,
         show_rownames = F)

pdf(paste("samplecorrel_heatmap.pdf",sep=''),height=18,width=22)
pheatmap
dev.off()

pcaData <- plotPCA(vsd, intgroup=c("Group", "Gender"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Group, shape=Gender)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_cowplot()
ggsave("fullPCA.pdf", width = 8, height = 6)

####Remove outliers and subset by sex
dds1 <- dds[,-c(7,45)] ##dds1 is the principle sample set

vsd1<-vst(dds1, blind=T)
sampleDists1 <- dist(t(assay(vsd1)))
sampleDistMatrix1 <- as.matrix(sampleDists1)
df1<-as.data.frame(colData(dds1)[,c("Group", "Gender")])
pheatmap1<-pheatmap(sampleDistMatrix1,
                   clustering_distance_rows=sampleDists1,
                   clustering_distance_cols=sampleDists1,
                   annotation_col=df1,
                   show_rownames=F)
pdf(paste("samplecorrel_removeout_heatmap.pdf",sep=''),height=18,width=22)
pheatmap
dev.off()

pcaData1 <- plotPCA(vsd1, intgroup=c("Group", "Gender"), returnData=TRUE)
percentVar1 <- round(100 * attr(pcaData1, "percentVar"))
ggplot(pcaData1, aes(PC1, PC2, color=Group, shape=Gender)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
  coord_fixed() +
  theme_cowplot()
ggsave("subsetPCA_fix.pdf", width = 8, height = 3)

pcaData2 <- plotPCA(vsd1, intgroup=c("Group", "age2"), returnData=TRUE)
percentVar2 <- round(100 * attr(pcaData2, "percentVar"))
ggplot(pcaData2, aes(PC1, PC2, color=Group, shape=age2)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  scale_shape_manual(values=c(0, 15))+
  coord_fixed() +
  theme_cowplot()
ggsave("subsetPCA_age.pdf", width = 8, height = 3)

#############Perform DESeq
dds1<-DESeq(dds1)
res1 <- results(dds1)
summary(res1)
##with adjusted p value <0.1; LFC>0 = 11, LFC<0 = 13, 41% of genes with low counts

res1.df<-as.data.frame(res1)
res1.df<-rownames_to_column(res1.df, var="ensembl")

genes.meta<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position", "chromosome_name", "entrezgene_id",
                               "external_gene_name", "gene_biotype"),
                  filters="ensembl_gene_id",
                  values=list(res1.df$ensembl),
                  mart=useMart("ensembl", dataset = "hsapiens_gene_ensembl"))

m<-match(res1.df$ensembl, genes.meta$ensembl_gene_id)
res1.df<-cbind(res1.df, genes.meta[m,])

write_csv(res1.df, "2022-12-04-MINI DESeq2 results_removed outliers.csv") ##annotation from getBM is better

###########################
##make heatmap of top 100 genes and volcano plot 

res1Ordered<-res1.df[order(abs(res1.df$stat), decreasing = T),]
res1Ordered<-res1Ordered %>% filter(!is.na(padj))
topgenes1<-res1Ordered[1:100,]
topgenes1[topgenes1 == ""]<-NA
topgenes1 <- topgenes1 %>% mutate(hgnc_symbol = coalesce(hgnc_symbol, ensembl))

vsd1.df<-as.data.frame(assay(vsd1))
vsd1.df<-rownames_to_column(vsd1.df, var = "orig.ensembl")
vsd1.df.top<-vsd1.df %>% filter(orig.ensembl %in% topgenes1$ensembl)
m<-match(vsd1.df.top$orig.ensembl, topgenes1$ensembl)
vsd1.df.top<-cbind(vsd1.df.top, topgenes1[m,])
rownames(vsd1.df.top)<-NULL
vsd1.df.top<-column_to_rownames(vsd1.df.top, var = "hgnc_symbol")
vsd1.df.top<-vsd1.df.top[2:80]

pheatmap3<- pheatmap(vsd1.df.top,
                    annotation_col=df1,
                    color = viridis(6),
                    main = "Top 100 DEG",
                    cluster_rows = T,
                    cluster_cols = T,
                    show_rownames = F,
                    show_colnames = F,
                    fontsize_col = 10, 
                    fontsize_row = 12,
                    scale = "row")
pdf(paste("genes_heatmap_nonames.pdf",sep=''),height=18,width=22)
pheatmap3
dev.off()

res1Ordered[res1Ordered==""]<-NA
res1Ordered<-res1Ordered %>% mutate(hgnc_symbol = coalesce(hgnc_symbol, ensembl))

EnhancedVolcano(res1Ordered, lab = res1Ordered$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                xlab=bquote(~Log[2]~ (frac("MD","Control"))),
                pCutoff = 0.1,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 10,
                labCol= "black", 
                boxedLabels = T,
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                subtitle="", 
                caption="",
                border="full",
                cutoffLineWidth=0.2,
                titleLabSize=10,
                ylim = c(0, 2),
                xlim = c(-2.5, 2.5))

ggsave2("MINI all volcano.pdf", device = "pdf", width = 12, height = 18)

###GSEA
res1.go<-res1.df[order(abs(res1.df$stat), decreasing = T),] %>% filter(!is.na(padj), !hgnc_symbol=="") %>% dplyr::select(hgnc_symbol, stat)
res1.go<-res1.go[!duplicated(res1.go$hgnc_symbol),]
geneList<-res1.go[,2]
names(geneList)=as.character(res1.go[,1])
geneList = sort(geneList, decreasing = T)
ego1<-gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "all", eps=0, pvalueCutoff = 0.05)
write_csv(ego1@result, "gseaGO_all.csv")

##format data for GOplot
GOres<-as.data.frame(ego1@result)
GOres<-GOres %>% dplyr::select("ONTOLOGY", "ID", "Description", "core_enrichment", "p.adjust")
names(GOres)<-c("Category", "ID", "Term", "Genes", "adj_pval")
rownames(GOres)<-NULL
GOres$Genes<-gsub("/", ", ", GOres$Genes)

GOdat<-res1.df %>% filter(!is.na(padj), !hgnc_symbol=="") %>% dplyr::select(hgnc_symbol, log2FoldChange, baseMean, stat, pvalue, padj)
names(GOdat)<-c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
rownames(GOdat)<-NULL

circ1<-circle_dat(GOres, GOdat)
reduced_circ1<-reduce_overlap(circ1, overlap = 0.75)
gobub1<-GOBubble(reduced_circ1, labels = 2, colour = c("#C2DF23", "#25858E", "#482173"))
ggsave("GOBubble_all.pdf", height = 15, width = 22) ##this drops the key; save manually to keep it

process<-c("structural constituent of ribosome",
           "mitochondrial protein-containing complex",
           "antigen binding",
           "immunoglobulin complex",
           "rRNA binding",
           "positive regulation of cell killing",
           "gated channel activity",
           "response to mechanical stimulus",
           "detection of chemical stimulus",
           "pattern recognition receptor activity",
           "G protein-coupled receptor activity",
           "positive regulation of interleukin-1 beta production")


circ_sub<-circ1 %>% filter(term %in% process)
genes2<-circ_sub %>% dplyr::select(genes, logFC)
names(genes2)<-c("ID", "logFC")

chord <- chord_dat(data = circ1, genes = genes2, process = process)

GOCluster(circ1, process, clust.by = 'logFC', term.col=viridis(12), lfc.col=c("red", "white", "dark blue"))
ggsave("GOCluster.pdf", width = 20, height = 12)

write_csv(ego1@result, "cluster profiler significant GO.csv")

########################################
##Figure 2

tm<-gmtPathways("~/Desktop/MINI-main/BTM_for_GSEA_20131008.gmt")

ranks1<-deframe(as.data.frame(res1.go))
fgseaRes1 <- fgsea(pathways=tm, stats=ranks1, maxSize=500)

collapsedPathways1 <- collapsePathways(fgseaRes1[order(pval)][padj < 0.01], 
                                      tm, ranks1)
mainPathways1 <- fgseaRes1[pathway %in% collapsedPathways1$mainPathways][
  order(-NES), pathway]
plotGseaTable(tm[mainPathways1], ranks1, fgseaRes1, gseaParam = 0.5) 

fgseaRes1$ID<-sub("^.*?\\(([M|S].*)\\)[^)]*$", "\\1", fgseaRes1$pathway)
fgseaRes1.2<-fgseaRes1 %>% rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ",")) %>% ungroup()
write_csv(fgseaRes1.2, "MD all fgsea.csv")
fgseaRes1_sig<-fgseaRes1[fgseaRes1$padj<0.05,]
fgseaRes1_sig<-fgseaRes1_sig %>% mutate(grp="MD")

########################################################################################################
##Import other datasets 

##Mende 2010 - GSE14882; GEO2R pipeline
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE14882", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000000000111111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("melas","control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
write_csv(tT, "./mende-melas-l2fc_2.csv")

melas.res<-tT %>% dplyr::select(Gene.symbol, t)
melas.res<-melas.res[!melas.res$Gene.symbol=="",]
melas.res<-melas.res[order(abs(melas.res$t), decreasing = T),] 
melas.res<-melas.res[!duplicated(melas.res$Gene.symbol),]

ranks2<-deframe(as.data.frame(melas.res))
fgsea_melas <- fgsea(pathways=tm, stats=ranks2, maxSize=500)
fgsea_melas$ID<-sub("^.*?\\(([M|S].*)\\)[^)]*$", "\\1", fgsea_melas$pathway)
fgsea_melas.2<-fgsea_melas %>% rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ",")) %>% ungroup()
write_csv(fgsea_melas.2, "melas all fgsea.csv")

m<-match(fgseaRes1_sig$ID, fgsea_melas$ID)
fgsea_melas_MD<-fgsea_melas[m,]
fgsea_melas_MD<-fgsea_melas_MD %>% mutate(grp="MELAS")

##########################################################################################
##West 2021
###mouse WT BMDM vs POLG BMDM - already in DESeq2 output (unstim vs unstim)
##direct import of these samples - avail on GEO - GSE171960

mouse.polg<-read_csv("WTBMDMunstim_vs_POLGBMDMunstim.genes.deseq.res.csv")
mouse.polg$GeneName<-sub("\\|.*", "\\1", mouse.polg$GeneName)
mouse.polg<-mouse.polg %>% mutate(GeneName = toupper(GeneName))
polg.res<-mouse.polg %>% dplyr::select(GeneName, stat)
polg.res<-na.omit(polg.res)
polg.res<-polg.res[order(abs(polg.res$stat), decreasing = T),]
polg.res<-polg.res[!duplicated(polg.res$GeneName),]
rankspolg<-deframe(as.data.frame(polg.res))

fgseaRespolg<-fgsea(pathways=tm, stats=rankspolg, maxSize=500)
fgseaRespolg$ID<-sub("^.*?\\(([M|S].*)\\)[^)]*$", "\\1", fgseaRespolg$pathway)
fgseaRespolg.2<-fgseaRespolg %>% rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ",")) %>% ungroup()
write_csv(fgseaRespolg.2, "WT-POLG BMDM fgsea res.csv")

m<-match(fgseaRes1_sig$ID, fgseaRespolg$ID)
fgseaPolg_MD<-fgseaRespolg[m,]
fgseaPolg_MD<-fgseaPolg_MD %>% mutate(grp="mouse Polg")

#################################################################################################################################
##West 2015 - GEO2R
##

# load series and platform data from GEO
gset <- getGEO("GSE63767", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "1100"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("case","control"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2)
tT <- topTable(fit2, adjust="fdr", sort.by = "B", number = "Inf")

write_csv(tT, "./west-tfam-l2fc.csv")

tfam<-tT %>% mutate(Gene.symbol = toupper(Gene.symbol))
tfam$Gene.symbol<-sub("\\|.*", "\\1", tfam$Gene.symbol)
tfam<-tfam %>% dplyr::select(Gene.symbol, t)

tfam <- tfam[!tfam$Gene.symbol=="",]
tfam<-tfam[order(abs(tfam$t), decreasing = T),]
tfam<-tfam[!duplicated(tfam$Gene.symbol),]

ranks.tfam<-deframe(as.data.frame(tfam))
fgseaRestfam<-fgsea(pathways=tm, stats=ranks.tfam, maxSize=500)
fgseaRestfam$ID<-sub("^.*?\\(([M|S].*)\\)[^)]*$", "\\1", fgseaRestfam$pathway)
fgseaRestfam.2<-fgseaRestfam %>% rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ",")) %>% ungroup()
write_csv(fgseaRestfam.2, "tfam MEF fgsea res.csv")

m<-match(fgseaRes1_sig$ID, fgseaRestfam$ID)
fgsea_tfam_MD<-fgseaRestfam[m,]
fgsea_tfam_MD<-fgsea_tfam_MD %>% mutate(grp="mouse Tfam")

#################################################################################################################################
allfgsea<-rbind(fgseaRes1_sig, fgsea_melas_MD, fgseaPolg_MD, fgsea_tfam_MD)

allfgsea <- allfgsea %>% mutate(Ann2 = case_when(
  (ID %in% c('M5.0','M95.0','M28','M71','M95.1','M200') ~ "Antigen Presentation"),
  (ID %in% c("M47.0", "M47.1", "M47.2", "M47.3", "M47.4","M54", "M69", "M83", "M9", "M58", "S2", "S8", "S9", "M156.0", "M156.1", "M123", "S3") ~ "B Cells and Plasma Cells"),
  (ID %in% c("M61.0", "M61.1", "M61.2", "M7.2", "M157", "S1") ~ "NK Cells"),
  (ID %in% c("M11.0", "M23", "M73", "M81", "M118.0","S4", "M4.15", "M118.1")~ "Monocytes"),
  (ID %in% c(" M4.9", "M7.0", "M7.1", "M7.4", "M14","M19", "M35.0", "M35.1", "M223", "M5.1", "M7.3", "M18", "M36","M44", "M52", "S0", "S6", "S7", "M4.11", "M4.5", "M4.9") ~ "T Cells"),
  (ID %in% c("M4.0", "M4.1", "M4.4", "M4.6","M4.7", "M4.10", "M6", "M8", "M22.0", "M22.1", "M31","M46", "M76", "M103", "M167", "M169", "M230", "M4.8", "M15","M37.3") ~ "Cell Cycle"),
  (ID %in% c("M40", "M43.0", "M43.1", "M53", "M64", "M67", "M86.0", "M86.1", "S5",  "S10", "M119", "M165", "M168") ~  "DC Cells"),
  (ID %in% c("M16", "M25", "M27.0", "M27.1", "M29","M37.0", "M115", "M146", "S11", "M33", "M50") ~ "Inflammatory TLR Chemokines"),
  (ID %in% c("M11.2", "M37.1", "M132", "M163") ~ "Neutrophils"),
  (ID %in% c("M75", "M127", "M158.0", "M158.1", "M13", "M68", "M111.0", "M111.1", "M150")~ "Interferon/Antiviral Sensing"),
  (ID %in% c("M187", "M238", "M231","M219","M216")~ "Mitochondria"),
  (ID %in% c("M87","M113","M37.2")~ "Golgi/ER"),
  (ID %in% c("M250", "M245", "M234", "M10.0", "M204.0", "M204.1")~"Transcription/Translation"),
  (ID %in% c("M3", "M4.2", "M59", "M130", "M4.3")~ "Signal Transduction"),
  
  (TRUE  ~ "0")
  
) )

write_csv(allfgsea, "BTM GSEA results all datasets_MD.csv")


p1<-ggplot(allfgsea, aes(y =reorder(pathway,-log10(padj)) , x = grp, size = -log10(padj), fill = NES)) + 
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
  guides(shape = guide_legend(override.aes = list(size = 8))) + 
  guides(color = guide_legend(override.aes = list(size = 8))) + 
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8)) + facet_grid(Ann2~.,scales="free", space="free_y")

pdf(paste("prelim NES plot 2.pdf"),height=13,width=7)
p1
dev.off()


##################################
##heatmaps of IFN plots

mouse.polg<-read_csv("WTBMDMunstim_vs_POLGBMDMunstim.genes.deseq.res.csv")
mouse.polg$GeneName<-sub("\\|.*", "\\1", mouse.polg$GeneName)
mouse.polg<-mouse.polg %>% mutate(GeneName = toupper(GeneName))

mouse.tfam<-read_csv("west-tfam-l2fc.csv")
mouse.tfam<-mouse.tfam %>% dplyr::select(Gene.symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
mouse.tfam<-mouse.tfam[order(abs(mouse.tfam$t), decreasing = T),]
mouse.tfam<-mouse.tfam[!duplicated(mouse.tfam$Gene.symbol),]
mouse.tfam$Gene.symbol<-sub("\\|.*", "\\1", mouse.tfam$Gene.symbol)
mouse.tfam<-mouse.tfam %>% mutate(Gene.symbol = toupper(Gene.symbol))

melas<-read_csv("mende-melas-l2fc_2.csv")
melas<-melas %>% dplyr::select(Gene.symbol, logFC, AveExpr, t, P.Value, adj.P.Val, B)
melas<-melas[order(abs(melas$t), decreasing = T),]
melas<-melas[!duplicated(melas$Gene.symbol),]

colnames(res1Ordered)
colnames(fgseaRes1_sig)

###clusters M75, M68

M75<-fgseaRes1_sig$leadingEdge[fgseaRes1_sig$ID=="M75"]
MD.m75<-res1Ordered %>% filter(hgnc_symbol %in% M75[[1]]) %>% dplyr::select(hgnc_symbol, log2FoldChange)
rownames(MD.m75)<-NULL
names(MD.m75)<-c("gene", "log2FC_MD")

melas.m75<-melas %>% filter(Gene.symbol %in% M75[[1]]) %>% dplyr::select(Gene.symbol, logFC)
rownames(melas.m75)<-NULL
names(melas.m75)<-c("gene", "log2FC_melas")
##CARD9 missing in melas 

polg.m75<-mouse.polg %>% filter(GeneName %in% M75[[1]]) %>% dplyr::select(GeneName, log2FoldChange)
rownames(polg.m75)<-NULL
names(polg.m75)<-c("gene", "log2FC_polg")

tfam.m75<-mouse.tfam %>% filter(Gene.symbol %in% M75[[1]]) %>% dplyr::select(Gene.symbol, logFC)
rownames(tfam.m75)<-NULL
names(tfam.m75)<-c("gene", "log2FC_tfam")


M75_l2fc<-merge(MD.m75, melas.m75, by.x = "gene") %>% merge(polg.m75, by.x = "gene") %>% merge(tfam.m75, by.x = "gene")
rownames(M75_l2fc)<-NULL
M75_l2fc<-column_to_rownames(M75_l2fc, var = "gene")
names(M75_l2fc)<-c("MD", "melas", "polg", "tfam")

breaks<-seq(-3,3, by=1) ###set consistent scale across plots 

p4<-pheatmap(M75_l2fc,
                     color = viridis(6),
                     main = "M75 cluster",
                     cluster_rows = T,
                     cluster_cols = F,
                     show_rownames = T,
                     show_colnames = T,
                     fontsize_col = 10, 
                     fontsize_row = 12,
                     breaks = breaks)

pdf(paste("M75_heatmap.pdf",sep=''),height=4,width=3)
p4
dev.off()



M68<-fgseaRes1_sig$leadingEdge[fgseaRes1_sig$ID=="M68"]
MD.M68<-res1Ordered %>% filter(hgnc_symbol %in% M68[[1]]) %>% dplyr::select(hgnc_symbol, log2FoldChange)
rownames(MD.M68)<-NULL
names(MD.M68)<-c("gene", "log2FC_MD")

melas.M68<-melas %>% filter(Gene.symbol %in% M68[[1]]) %>% dplyr::select(Gene.symbol, logFC)
rownames(melas.M68)<-NULL
names(melas.M68)<-c("gene", "log2FC_melas")

polg.M68<-mouse.polg %>% filter(GeneName %in% M68[[1]]) %>% dplyr::select(GeneName, log2FoldChange)
rownames(polg.M68)<-NULL
names(polg.M68)<-c("gene", "log2FC_polg")

tfam.M68<-mouse.tfam %>% filter(Gene.symbol %in% M68[[1]]) %>% dplyr::select(Gene.symbol, logFC)
rownames(tfam.M68)<-NULL
names(tfam.M68)<-c("gene", "log2FC_tfam")
 
M68_l2fc<-merge(MD.M68, melas.M68, by.x = "gene") %>% merge(polg.M68, by.x = "gene") %>% merge(tfam.M68, by.x = "gene")
rownames(M68_l2fc)<-NULL
M68_l2fc<-column_to_rownames(M68_l2fc, var = "gene")
names(M68_l2fc)<-c("MD", "melas", "polg", "tfam")

p5<-pheatmap(M68_l2fc,
             color = viridis(6),
             main = "M68 cluster",
             cluster_rows = T,
             cluster_cols = F,
             show_rownames = T,
             show_colnames = T,
             fontsize_col = 10, 
             fontsize_row = 12,
             breaks = breaks)

pdf(paste("M68_heatmap.pdf",sep=''),height=4,width=3)
p5
dev.off()






