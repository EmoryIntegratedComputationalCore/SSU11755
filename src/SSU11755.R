pacman::p_load("here", "readr", "janitor", "assertr", "tidyverse",
               "DESeq2", "ggplot2", "EnhancedVolcano", "apeglm", "vsn",
               "pheatmap", "org.Mm.eg.db", "AnnotationDbi",
               "pathview", "gage", "gageData", "dplyr")
data(kegg.sets.mm)
theme_set(theme_classic())

files <- list(
  counts = here("SSU11755/Repository/input/countdata.csv"),
  sampledata = here("SSU11755/Repository/input/sampledata.csv"),
  results_res1 = here("SSU11755/Repository/output/SSU11755_6vsNaive_DEGanalysis.csv"),
  results_res2 = here("SSU11755/Repository/output/SSU11755_9vsNaive_DEGanalysis.csv"),
  results_res3 = here("SSU11755/Repository/output/SSU11755_9vs6_DEGanalysis.csv"),
  graph_pca = here("SSU11755/Repository/output/graphs/SSU11755_pca_group.png"),
  graph_ma1 = here("SSU11755/Repository/output/graphs/SSU11755_ma_6vsNaive.png"),
  graph_ma2 = here("SSU11755/Repository/output/graphs/SSU11755_ma_9vsNaive.png"),
  graph_ma3 = here("SSU11755/Repository/output/graphs/SSU11755_ma_9vs6.png"),
  graph_volc1 = here("SSU11755/Repository/output/graphs/SSU11755_volcano_6vsNaive.png"),
  graph_volc2 = here("SSU11755/Repository/output/graphs/SSU11755_volcano_9vsNaive.png"),
  graph_volc3 = here("SSU11755/Repository/output/graphs/SSU11755_volcano_9vs6.png"),
  graph_heat1 = here("SSU11755/Repository/output/graphs/SSU11755_heatmapNaiveref.png"),
  graph_heat2 = here("SSU11755/Repository/output/graphs/SSU11755_heatmapNaiveref2.png"),
  graph_heat3 = here("SSU11755/Repository/output/graphs/SSU11755_heatmap6ref.png")
  )

#load count and sample data for all comparisons as a glm
countdata <- as.matrix(read.csv(files$counts, row.names=1))

#sample data
sampledata <- read.csv(files$sampledata, row.names=1) %>%
  janitor::clean_names()

stopifnot(rownames(sampledata) %in% colnames(countdata))

dds1 <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampledata,
                              design = ~batch * dose_group)
#set the reference level
#comparing B to A and C to A
dds1$dose_group <- relevel(dds1$dose_group, ref = "Naïve")
dds1 <- DESeq(dds1)

#6 to Naive
res1 <- results(dds1, alpha=0.01, name="dose_group_1x10.6.HAd_NP_vs_Naïve")
summary(res1)

#9 to naive
res2 <- results(dds1, alpha=0.01, name="dose_group_1x10.9.HAd_NP_vs_Naïve")
summary(res2)

#create data frames of results of each comparison with gene symbol in a column called  ID
res1_df <- as.data.frame(res1) %>%
  mutate(ID = as.factor((row.names(res1))))
stopifnot(nrow(res1_df) == 26485 & ncol(res1_df) == 7)

res2_df <- as.data.frame(res2) %>%
  mutate(ID = as.factor((row.names(res2))))
stopifnot(nrow(res2_df) == 26485 & ncol(res2_df) == 7)

#export
res1_df %>%
  write_delim(files$results_res1, delim=",")

res2_df %>%
  write_delim(files$results_res2, delim=",")

####PCA
vsd1 <- vst(dds1, blind=FALSE)

(pca <- plotPCA(vsd1, intgroup=c("dose_group")))

#export
ggsave(files$graph_pca, plot=pca, dpi=700)

####MA #FIXME: change to most recent coef names
#compare 6 vs Naive
resLFC1 <- lfcShrink(dds1, coef="dose_group_1x10.6.HAd_NP_vs_Naïve", type="apeglm")

#compare 9 vs Naive
resLFC2 <- lfcShrink(dds1, coef="dose_group_1x10.9.HAd_NP_vs_Naïve", type="apeglm")

#plot 1
plotMA(resLFC1, ylim=c(-10, 10), main="DE genes between 1x10^6 and Naïve")

maplot1 <- png(files$graph_ma1, width=450, height=450)
plotMA(resLFC1, ylim=c(-10, 10), main="DE genes between 1x10^6 and Naïve")
dev.off()

#plot 2
plotMA(resLFC2, ylim=c(-10, 10), main="DE genes between 1x10^9 and Naïve")

maplot2 <- png(files$graph_ma2, width=450, height=450)
plotMA(resLFC2, ylim=c(-10, 10), main="DE genes between 1x10^9 and Naïve")
dev.off()

###Volcano
#res1
volc1 <- EnhancedVolcano(res1,
                         lab = rownames(res1),
                         x = "log2FoldChange",
                         y = "padj",
                         xlim = c(-10, 10),
                         title= NULL,
                         subtitle= NULL,
                         pLabellingCutoff = 0.01,
                         pCutoff = 0.01,
                         legendPosition = "bottom",
                         legend=c("NS", "Log2 fold-change", "adj P-value",
                                  "adj P-value & Log2 fold-change"))

volc1 <- png(files$graph_volc1, width=700, height=600)
EnhancedVolcano(res1,
                lab = rownames(res1),
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-10, 10),
                title= NULL,
                subtitle= NULL,
                pLabellingCutoff = 0.01,
                pCutoff = 0.01,
                legendPosition = "bottom",
                legend=c("NS", "Log2 fold-change", "adj P-value",
                         "adj P-value & Log2 fold-change"))
dev.off()

#res2
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-10, 10),
                title= NULL,
                subtitle= NULL,
                pLabellingCutoff = 0.01,
                pCutoff = 0.01,
                legendPosition = "bottom",
                legend=c("NS", "Log2 fold-change", "adj P-value",
                         "adj P-value & Log2 fold-change"))

volc2 <- png(files$graph_volc2, width=700, height=600)
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-10, 10),
                title=NULL,
                subtitle=NULL,
                pLabellingCutoff = 0.01,
                pCutoff = 0.01,
                legendPosition = "bottom",
                legend=c("NS", "Log2 fold-change", "adj P-value",
                         "adj P-value & Log2 fold-change"))
dev.off()

###Heatmaps
#one for each comparison using vsd transformed data created for PCA plots

#1, 6 vs Naive

#remove samples from batch 1
vsd1_filt <- as.data.frame(assays(vsd1))
vsd1_filt[c(1:5)] <- NULL

#remove samples in 9 group
vsd1_filt[2] <- NULL
vsd1_filt[4] <- NULL

#reorder by dose group from left to right
vsd1_filt <- vsd1_filt[c(2, 4, 1, 3)]

#select only top 20 DE genes
select1 <- order(res1_df$padj, decreasing = FALSE)[1:20]

#select vars of interest
heat1 <- as.data.frame(colData(dds1)["dose_group"])

#change order in which columns appear
callback = function(hc,mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram (hc), wts = sv)
  as.hclust(dend)
}

#specify naive as green, 6 is pink
#set color order
colororder1 = list(
  dose_group = c("Naïve"="#ccebc5", "1x10^6 HAd_NP"="#fbb4ae"))

#plot
ht1 <- pheatmap(vsd1_filt[select1, ],
                cluster_rows=FALSE,
                show_rownames = TRUE,
                cluster_cols = FALSE,
                annotation_col = heat1,
                clustering_callback = callback,
                annotation_colors = colororder1,
                cellwidth = 30,
                cellheight = 8,
                fontsize = 9)

#export
ggsave(files$graph_heat1, plot=ht1, dpi=600)

#2, 9 vs Naive
#remove samples in batch 1 group
vsd2_filt <- as.data.frame(assays(vsd1))
vsd2_filt[c(1:5)] <- NULL

#remove samples in 6 group
vsd2_filt[1] <- NULL
vsd2_filt[3] <- NULL

#reorder by dose group from left to right
vsd2_filt <- vsd2_filt[c(2, 4, 1, 3)]

#top 20 DE genes in comparison
select2 <- order(res2_df$padj, decreasing = FALSE)[1:20]

#choose variable of interest
heat2 <- as.data.frame(colData(dds1)["dose_group"])

#specify naive as green, 9 is blue
#set color order
colororder2 = list(
  dose_group = c("Naïve"="#ccebc5", "1x10^9 HAd_NP"="#b3cde4"))

#plot
ht2 <- pheatmap(vsd2_filt[select2, ], 
                cluster_rows=FALSE, 
                show_rownames = TRUE,
                cluster_cols = FALSE, 
                annotation_col = heat2, 
                clustering_callback = callback,
                annotation_colors = colororder2,
                cellwidth = 30,
                cellheight = 8,
                fontsize = 9)

#export
ggsave(files$graph_heat2, plot=ht2, dpi=600)

#comparing 9 to 6
dds1$dose_group <- relevel(dds1$dose_group, ref = "1x10^6 HAd_NP")
dds1 <- DESeq(dds1)

#9 to 6 
res3 <- results(dds1, alpha=0.01, name="dose_group_1x10.9.HAd_NP_vs_1x10.6.HAd_NP")
summary(res3)

#create data frames of results of each comparison with gene symbol in a column called  ID
res3_df <- as.data.frame(res3) %>%
  mutate(ID = as.factor((row.names(res3))))
stopifnot(nrow(res3_df) == 26485 & ncol(res3_df) == 7)

res3_df %>%
  write_delim(files$results_res3, delim=",")

####MA
#compare 9 vs 6
resLFC3 <- lfcShrink(dds1, coef="dose_group_1x10.9.HAd_NP_vs_1x10.6.HAd_NP", type="apeglm")

#plot 3
plotMA(resLFC3, ylim=c(-5, 5), main="DE genes between 1x10^9 and 1x10^6")

maplot3 <- png(files$graph_ma3, width=450, height=450)
plotMA(resLFC3, ylim=c(-5, 5), main="DE genes between 1x10^9 and 1x10^6")
dev.off()

#res3
EnhancedVolcano(res3,
                lab = rownames(res3),
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-10, 10),
                title=NULL,
                subtitle=NULL,
                pLabellingCutoff = 0.01,
                pCutoff = 0.01,
                legendPosition = "bottom",
                legend=c("NS", "Log2 fold-change", "adj P-value",
                         "adj P-value & Log2 fold-change"))

volc3 <- png(files$graph_volc3, width=700, height=600)
EnhancedVolcano(res3,
                lab = rownames(res3),
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-10, 10),
                title=NULL,
                subtitle=NULL,
                pLabellingCutoff = 0.01,
                pCutoff = 0.01,
                legendPosition = "bottom",
                legend=c("NS", "Log2 fold-change", "adj P-value",
                         "adj P-value & Log2 fold-change"))
dev.off()

###Heatmaps
#3, 9 vs 6

#VST transformation
vsd3 <- vst(dds1, blind=FALSE)

#remove samples from batch 1
vsd3_filt <- as.data.frame(assays(vsd3))
vsd3_filt[c(1:5)] <-NULL

#remove samples in naive group
vsd3_filt[3] <- NULL
vsd3_filt[6] <- NULL

#reorder by dose group from left to right
vsd3_filt <- vsd3_filt[c(2, 4, 3, 1)]

#top 20 DE genes in comparison
select3 <- order(res3_df$padj, decreasing = FALSE)[1:20]

#choose variable of interest
heat3 <- as.data.frame(colData(dds1)["dose_group"])

#specify naive as green
#set color order
colororder3 = list(
  dose_group = c("1x10^6 HAd_NP"="#fbb4ae", "1x10^9 HAd_NP"="#b3cde4"))

ht3<-pheatmap(vsd3_filt[select3,], 
              cluster_rows=FALSE,
              show_rownames = TRUE,
              cluster_cols = FALSE,
              annotation_col = heat3,
              clustering_callback = callback,
              annotation_colors = colororder3,
              cellwidth = 30,
              cellheight = 8,
              fontsize = 9)

#export
ggsave(files$graph_heat3,plot=ht3,dpi=600)

#############################
#pathway analysis
#############################
#get entrez ids and genenames
res1$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res1),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

res1$name <- mapIds(org.Mm.eg.db,
                     keys=row.names(res1),
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")

head(res1, 10)

#res2
res2$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res2),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

res2$name <- mapIds(org.Mm.eg.db,
                     keys=row.names(res2),
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")

head(res2, 10)

#res3
res3$entrez <- mapIds(org.Mm.eg.db,
                     keys=row.names(res3),
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

res3$name <- mapIds(org.Mm.eg.db,
                     keys=row.names(res3),
                     column="GENENAME",
                     keytype="SYMBOL",
                     multiVals="first")

head(res3, 10)

#KEGG pathway analysis

#get a matrix of fold changes and entrez ids
foldchanges1 <- res1$log2FoldChange
names(foldchanges1) <- res1$entrez

#res2
foldchanges2 <- res2$log2FoldChange
names(foldchanges2) <- res2$entrez

#res3
foldchanges3 <- res3$log2FoldChange
names(foldchanges3) <- res3$entrez

#run KEGG pathway analysis
keggres1 <- gage(foldchanges1, gsets= kegg.sets.mm, same.dir=TRUE)

# Look at both up (greater), down (less), and statistics
lapply(keggres1, head, 10)

#no pathways at FDR <0.1, complement and coagulation at 0.2

#keggres2
#no DE pathways

#keggres3
keggres3 <- gage(foldchanges3, gsets= kegg.sets.mm, same.dir=TRUE)
# Look at both up (greater), down (less), and statistics
lapply(keggres3, head, 10)
#same DE pathway at FDR 0.12