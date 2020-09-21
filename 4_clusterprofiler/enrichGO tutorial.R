
######################
a <- bitr(x , fromType = "SYMBOL",
                toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)

GO1 <- enrichGO(gene          = a$ENTREZID, 
                  universe      = uni$ENTREZID, 
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP", #CC #MF
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
GO1 <- setReadable(GO1, OrgDb = org.Hs.eg.db)
summary(GO1)

##########Filter
g1<- dropGO(nego1, level = NULL, term = NULL) ###remove specific go 
g2<- gofilter(nego1, level = 4)
g3<- simplify(nego1, cutoff = 0.7, by = "p.adjust", select_fun=min,
              measure = "Wang", semData = NULL)

###########Graph
barplot(GO1)

p1 <- dotplot(GO1, showCategory=10)
p1
p2 <- dotplot(GO1, showCategory=10)
p2
plot_grid(p1, p2, ncol=2)

cnetplot(GO1, foldChange = geneList)
cnetplot(GO1, foldChange = geneList, categorySize="pvalue")
cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)

heatplot(edox)
heatplot(edox, foldChange=geneList)

upsetplot(edo)