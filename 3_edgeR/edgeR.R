#edgeR

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")
library("edgeR")
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)

##START
load("tcga.RData")
pheno <- as.data.frame(colData(dataPrep1))
htseqcounts <- assay(dataPrep1)
rm( query, queryDown, dataSmTM,dataSmTP, DataDirectory, FileNameData,
   rownames, samplesDown, gene.location)

###########
htseqcounts<- t(htseqcounts)
htseqcounts <- as.data.frame(htseqcounts)

Percentile_00  = quantile(htseqcounts$ENSG00000108064,0)
Percentile_25  = quantile(htseqcounts$ENSG00000108064, 0.25)
Percentile_75  = quantile(htseqcounts$ENSG00000108064, 0.75)
Percentile_100 = quantile(htseqcounts$ENSG00000108064,1)


RB = rbind(Percentile_00, Percentile_25, Percentile_75, Percentile_100)

dimnames(RB)[[2]] = "Value"

RB

htseqcounts$Group[htseqcounts$ENSG00000108064 >= Percentile_00 & htseqcounts$ENSG00000108064 <  Percentile_25]  = "low"
htseqcounts$Group[htseqcounts$ENSG00000108064 >= Percentile_25 & htseqcounts$ENSG00000108064 <  Percentile_75]  = "middle"
htseqcounts$Group[htseqcounts$ENSG00000108064 >= Percentile_75 & htseqcounts$ENSG00000108064 <= Percentile_100] = "high"

htseqcounts$Group = factor(htseqcounts$Group,
                           levels=c("high", "middle", "low"))

#####################
group <- htseqcounts$Group
htseqcounts <- assay(dataPrep1)
htseqcounts <- as.data.frame(htseqcounts)
#
y <- DGEList(counts=htseqcounts, group = group)
levels(y$samples$group)
head(y$samples)
#
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
#
y <- calcNormFactors(y)
head(y$samples)
#classic
y <- estimateDisp(y)
et <- exactTest(y, pair=c("high", "low"))
topTags(et) #adjusted p-value
dgeresult <- et$table

############
dgeresult <- cbind(rownames(dgeresult), dgeresult)
fil <- dgeresult %>% filter(PValue < 0.01)
filneg <- fil %>% filter(logFC < -1.5)
filpos <- fil %>% filter(logFC > 1.5)
dgefiltered <- rbind(filneg, filpos)
rm(data, design, et, fil,  htseqcounts, pheno, RB, y ,
   group, keep, Percentile_00, Percentile_100, Percentile_25, Percentile_75)
