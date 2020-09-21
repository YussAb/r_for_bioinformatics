library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(stringr)
library(DOSE)
library(ReactomePA)
library(biomaRt)
library(tidyverse)

#https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html

uni <- read.csv("nature.csv")
universe <- as.character(uni$Entrez.Gene.ID)

###################
met <- as.data.frame(uni$Gene.Symbol)

DEfiltered1 <- filneg
Defneg1 = DEfiltered1[,1]
Defneg1 <- bitr(Defneg1 , fromType = "ENSEMBL",
             toType = c("ENTREZID", "SYMBOL"),
             OrgDb = org.Hs.eg.db)
colnames(met) <- c("SYMBOL")
neg<- merge(Defneg1, met)

#2
DEfiltered1 <- filpos
Defpos1 = DEfiltered1[,1]
Defpos1 <- bitr(Defpos1 , fromType = "ENSEMBL",
             toType = c("ENTREZID", "SYMBOL"),
             OrgDb = org.Hs.eg.db)
pos <- merge(Defpos1, met)

###Enrichment Pathway 

ego <- rbind(neg, pos)

################
library(ReactomePA)

x <- enrichPathway(gene=neg$ENTREZID,pvalueCutoff=0.05, readable=T)
y <- enrichPathway(gene=pos$ENTREZID,pvalueCutoff=0.05, readable=T)
z <- enrichPathway(gene=ego$ENTREZID,pvalueCutoff=0.05, readable=T)

emapplot(x)

emapplot(y)

emapplot(z)
