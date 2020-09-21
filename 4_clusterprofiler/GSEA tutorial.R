###GSEA
#4 GENE ENRICHMENT ANALYSIS
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(DOSE)
library(enrichplot)

d = read.csv(your_csv_file)
## assume 1st column is ID
## 2nd column is FC
## feature 1: numeric vector
geneList = d[,2]
## feature 2: named vector
names(geneList) = as.character(d[,1])
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

##
gsecc1 <- gseGO(geneList=Def1, ont="CC", OrgDb=org.Mm.eg.db, verbose=F)
head(summary(gsecc1))
a<- gsecc1@result

##Graph
gseaplot2(gsecc1, geneSetID="GO:0044233", pvalue_table = T, ES_geom = "line" )
gseaplot(gsecc1, geneSetID="GO:0006979", title= "Oxidative Stress")

dotplot(gsecc1, showCategory=50)
cnetplot(gsecc1, foldChange= cor)
heatplot(gsecc1, foldChange= cor)
emapplot(gsecc1)
upsetplot(gsecc1)
ridgeplot(gsecc1)
