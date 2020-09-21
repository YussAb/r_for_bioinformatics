https://yulab-smu.github.io/clusterProfiler-book/chapter7.html
#############KEGG

kk <- enrichKEGG(gene         = Defpos1$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk, OrgDb = org.Hs.eg.db)
head(kk)
dotplot(kk)


###############
kk2 <- gseKEGG(geneList     = cor,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
dotplot(kk2)
