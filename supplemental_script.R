####supplemental figures
##requires packages and some objects from main text script

ddsm <- dds1[, dds1$Gender == "M"]
ddsf <- dds1[, dds1$Gender == "F"]


#######male samples 

vsdm<-vst(ddsm, blind=T)
sampleDistsm <- dist(t(assay(vsdm)))
sampleDistMatrixm <- as.matrix(sampleDistsm)
dfm<-as.data.frame(colData(ddsm)[,c("Group", "Gender")])
pheatmapm<-pheatmap(sampleDistMatrixm,
                    clustering_distance_rows=sampleDistsm,
                    clustering_distance_cols=sampleDistsm,
                    annotation_col=dfm,
                    show_rownames=F,
                    show_colnames=F)
pdf(paste("samplecorrel_male_heatmap.pdf",sep=''),height=10,width=12)
pheatmapm
dev.off()

pcaDatam <- plotPCA(vsdm, intgroup=c("Group"), returnData=TRUE)
percentVarm <- round(100 * attr(pcaDatam, "percentVar"))
ggplot(pcaDatam, aes(PC1, PC2, color=Group)) +
  geom_point(size=3, shape = 17) +
  xlab(paste0("PC1: ",percentVarm[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarm[2],"% variance")) + 
  coord_fixed() +
  theme_cowplot()
ggsave("male PCA.pdf", width = 8, height = 6)

ddsm<-DESeq(ddsm)
resm <- results(ddsm)
summary(resm)
##with adjusted p value <0.1; LFC>0 = 6, LFC<0 = 12, 62% of genes with low counts

resm.df<-as.data.frame(resm)
resm.df<-rownames_to_column(resm.df, var="ensembl")

m<-match(resm.df$ensembl, genes.meta$ensembl_gene_id)
resm.df<-cbind(resm.df, genes.meta[m,])
write_csv(resm.df, "2022-12-04-MINI DESeq2 results_male.csv")

resmOrdered<-resm.df[order(abs(resm.df$stat), decreasing = T),]
resmOrdered<-resmOrdered %>% filter(!is.na(padj))
topgenesm<-resmOrdered[1:100,]
topgenesm[topgenesm == ""]<-NA
topgenesm <- topgenesm %>% mutate(hgnc_symbol = coalesce(hgnc_symbol, ensembl))

vsdm.df<-as.data.frame(assay(vsdm))
vsdm.df<-rownames_to_column(vsdm.df, var = "orig.ensembl")
vsdm.df.top<-vsdm.df %>% filter(orig.ensembl %in% topgenesm$ensembl)
m<-match(vsdm.df.top$orig.ensembl, topgenesm$ensembl)
vsdm.df.top<-cbind(vsdm.df.top, topgenesm[m,])
rownames(vsdm.df.top)<-NULL
vsdm.df.top<-column_to_rownames(vsdm.df.top, var = "hgnc_symbol")
vsdm.df.top<-vsdm.df.top[,2:39]

pheatmap_m<- pheatmap(vsdm.df.top,
                     annotation_col=dfm,
                     color = viridis(6),
                     main = "Top 100 DEG",
                     cluster_rows = T,
                     cluster_cols = T,
                     show_rownames = T,
                     show_colnames = F,
                     fontsize_col = 10, 
                     fontsize_row = 12,
                     scale = "row")

pdf(paste("genes_heatmap_male.pdf",sep=''),height=18,width=22)
pheatmap_m
dev.off()

resmOrdered[resmOrdered==""]<-NA
resmOrdered<-resmOrdered %>% mutate(hgnc_symbol = coalesce(hgnc_symbol, ensembl))

EnhancedVolcano(resmOrdered, lab = resmOrdered$hgnc_symbol,
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

ggsave2("MINI male volcano.pdf", device = "pdf", width = 12, height = 18)


EnhancedVolcano(resmOrdered, lab = resmOrdered$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                xlab=bquote(~Log[2]~ (frac("MD","Control"))),
                pCutoff = 0.1,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 4,
                labCol= "black", 
                boxedLabels = T,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                subtitle="", 
                caption="",
                border="full",
                cutoffLineWidth=0.2,
                titleLabSize=10,
                ylim = c(0, 3.5),
                xlim = c(-2.5, 2.5))

ggsave2("MINI male volcano.pdf", device = "pdf", width = 12, height = 18)

###GSEA
resm.go<-resm.df[order(abs(resm.df$stat), decreasing = T),] %>% filter(!is.na(padj), !hgnc_symbol=="") %>% dplyr::select(hgnc_symbol, stat)
resm.go<-resm.go[!duplicated(resm.go$hgnc_symbol),]
geneList<-resm.go[,2]
names(geneList)=as.character(resm.go[,1])
geneList = sort(geneList, decreasing = T)
egom<-gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "all", eps=0, pvalueCutoff = 0.05)
write_csv(egom@result, "gseaGO male.csv")

##format data for GOplot
GOresm<-as.data.frame(egom@result)
GOresm<-GOresm %>% dplyr::select("ONTOLOGY", "ID", "Description", "core_enrichment", "p.adjust")
names(GOresm)<-c("Category", "ID", "Term", "Genes", "adj_pval")
rownames(GOresm)<-NULL
GOresm$Genes<-gsub("/", ", ", GOresm$Genes)

GOdatm<-resm.df %>% filter(!is.na(padj), !hgnc_symbol=="") %>% dplyr::select(hgnc_symbol, log2FoldChange, baseMean, stat, pvalue, padj)
names(GOdatm)<-c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
rownames(GOdatm)<-NULL

circm<-circle_dat(GOresm, GOdatm)
reduced_circm<-reduce_overlap(circm, overlap = 0.75)
gobubm<-GOBubble(reduced_circm, labels = 2, colour = c("#C2DF23", "#25858E", "#482173"))
pdf(paste("GObub_male.pdf",sep=''),height=3,width=8)

ggsave("GOBub_male.pdf", height = 15, width = 22)

#############
##female samples

vsdf<-vst(ddsf, blind=T)
sampleDistsf <- dist(t(assay(vsdf)))
sampleDistMatrixf <- as.matrix(sampleDistsf)
dff<-as.data.frame(colData(ddsf)[,c("Group", "Gender")])
pheatmapf<-pheatmap(sampleDistMatrixf,
                    clustering_distance_rows=sampleDistsf,
                    clustering_distance_cols=sampleDistsf,
                    annotation_col=dff,
                    show_rownames=F,
                    show_colnames=F)
pdf(paste("samplecorrel_female_heatmap.pdf",sep=''),height=10,width=12)
pheatmapf
dev.off()

pcaDataf <- plotPCA(vsdf, intgroup=c("Group"), returnData=TRUE)
percentVarf <- round(100 * attr(pcaDataf, "percentVar"))
ggplot(pcaDataf, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVarf[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarf[2],"% variance")) + 
  coord_fixed() +
  theme_cowplot()
ggsave("female PCA.pdf", width = 8, height = 6)

ddsf<-DESeq(ddsf)
resf <- results(ddsf)
summary(resf)
##with adjusted p value <0.1; LFC>0 = 0, LFC<0 = 0

resf.df<-as.data.frame(resf)
resf.df<-rownames_to_column(resf.df, var="ensembl")

m<-match(resf.df$ensembl, genes.meta$ensembl_gene_id)
resf.df<-cbind(resf.df, genes.meta[m,])
write_csv(resf.df, "2022-12-04-MINI DESeq2 results_female.csv")

resfOrdered<-resf.df[order(abs(resf.df$stat), decreasing = T),]
resfOrdered<-resfOrdered %>% filter(!is.na(padj))
topgenesf<-resfOrdered[1:100,]
topgenesf[topgenesf == ""]<-NA
topgenesf <- topgenesf %>% mutate(hgnc_symbol = coalesce(hgnc_symbol, ensembl))

vsdf.df<-as.data.frame(assay(vsdf))
vsdf.df<-rownames_to_column(vsdf.df, var = "orig.ensembl")
vsdf.df.top<-vsdf.df %>% filter(orig.ensembl %in% topgenesf$ensembl)
m<-match(vsdf.df.top$orig.ensembl, topgenesf$ensembl)
vsdf.df.top<-cbind(vsdf.df.top, topgenesf[m,])
rownames(vsdf.df.top)<-NULL
vsdf.df.top<-column_to_rownames(vsdf.df.top, var = "hgnc_symbol")
vsdf.df.top<-vsdf.df.top[,2:42]

pheatmap_f<- pheatmap(vsdf.df.top,
                      annotation_col=dff,
                      color = viridis(6),
                      main = "Top 100 DEG",
                      cluster_rows = T,
                      cluster_cols = T,
                      show_rownames = T,
                      show_colnames = F,
                      fontsize_col = 10, 
                      fontsize_row = 12,
                      scale = "row")

pdf(paste("genes_heatmap_female.pdf",sep=''),height=18,width=22)
pheatmap_f
dev.off()

resfOrdered[resfOrdered==""]<-NA
resfOrdered<-resfOrdered %>% mutate(hgnc_symbol = coalesce(hgnc_symbol, ensembl))

EnhancedVolcano(resfOrdered, lab = resfOrdered$hgnc_symbol,
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

ggsave2("MINI female volcano.pdf", device = "pdf", width = 12, height = 18)

EnhancedVolcano(resfOrdered, lab = resfOrdered$hgnc_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                xlab=bquote(~Log[2]~ (frac("MD","Control"))),
                pCutoff = 0.1,
                FCcutoff = 0.5,
                pointSize = 4,
                labSize = 4,
                labCol= "black", 
                boxedLabels = T,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30',
                subtitle="", 
                caption="",
                border="full",
                cutoffLineWidth=0.2,
                titleLabSize=10,
                ylim = c(0, 3.5),
                xlim = c(-2.5, 2.5))

ggsave2("MINI female volcano.pdf", device = "pdf", width = 12, height = 18)

###GSEA
resf.go<-resf.df[order(abs(resf.df$stat), decreasing = T),] %>% filter(!is.na(padj), !hgnc_symbol=="") %>% dplyr::select(hgnc_symbol, stat)
resf.go<-resf.go[!duplicated(resf.go$hgnc_symbol),]
geneList<-resf.go[,2]
names(geneList)=as.character(resf.go[,1])
geneList = sort(geneList, decreasing = T)
egof<-gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "all", eps=0, pvalueCutoff = 0.05)
write_csv(egof@result, "gseaGO female.csv")

##format data for GOplot
GOresf<-as.data.frame(egof@result)
GOresf<-GOresf %>% dplyr::select("ONTOLOGY", "ID", "Description", "core_enrichment", "p.adjust")
names(GOresf)<-c("Category", "ID", "Term", "Genes", "adj_pval")
rownames(GOresf)<-NULL
GOresf$Genes<-gsub("/", ", ", GOresf$Genes)

GOdatf<-resf.df %>% filter(!is.na(padj), !hgnc_symbol=="") %>% dplyr::select(hgnc_symbol, log2FoldChange, baseMean, stat, pvalue, padj)
names(GOdatf)<-c("ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
rownames(GOdatf)<-NULL

circf<-circle_dat(GOresf, GOdatf)
reduced_circf<-reduce_overlap(circf, overlap = 0.75)
gobubf<-GOBubble(reduced_circf, labels = 2, colour = c("#C2DF23", "#25858E", "#482173"))
ggsave("GOBubble_female.pdf", height = 15, width = 22)









