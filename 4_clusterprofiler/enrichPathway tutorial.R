#########Enrich Pathway

x <- enrichPathway(gene=filneg$ENTREZID, pvalueCutoff=0.05, readable=T)
emapplot(x)