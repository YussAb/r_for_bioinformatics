library(Biobase)
library(GEOquery)
library(limma)


#data_from_article
load('GSE72267_data.RData')


#function
gse_pt_class <- function(x){
  library(Biobase)
  library(GEOquery)
  library(limma)
  a <- getGEO(x, GSEMatrix =TRUE, AnnotGPL=TRUE)
  exp_data <- as.data.frame(exprs(a[[1]]))
  pheno_data <- a[[paste(x,'_series_matrix.txt.gz', sep='')]]@phenoData@data
  patient_class <- cbind.data.frame(geo_accession = pheno_data$geo_accession, title=pheno_data$title)
  exp_data <- as.data.frame(t(exp_data)) #block
  exp_data <- cbind.data.frame(geo_accession = rownames(exp_data),exp_data)
  y <-merge(exp_data,patient_class)
  row.names(y) <- y$title
  y$title <-NULL
  output <- as.data.frame(t(y))
  return(output)
}

gse72267 <- gse_pt_class('GSE72267') #call function

gse72267 <- gse72267[-1,] #remove control 

############http://biolearnr.blogspot.com/2017/05/bfx-clinic-getting-up-to-date.html#######
GEOset <- GSE72267[[1]]
featureData(GEOset)

BiocManager::install('hgu133a2.db')
library('hgu133a2.db')
library(dplyr)

keys(hgu133a2.db)
(rownames(GEOset)) %in% keys(hgu133a2.db) %>%  summary
columns(hgu133a2.db)

###########
ae.annots <- AnnotationDbi::select(
  x       = hgu133a2.db,
  keys    = rownames(GEOset),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
dim(ae.annots)

# Redefinition of ae.annots
ae.annots <- AnnotationDbi::select(
  x       = hgu133a2.db,
  keys    = rownames(GEOset),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
) %>%
  group_by(PROBEID) %>%
  summarise_each(funs(collapser)) %>%
  ungroup



#############https://www.biostars.org/p/332461/#############
require(GEOquery)
require(Biobase)
#gset <- GSE72267
#if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
#gset <- gset[[idx]]
#annotlookup
load("~/UNITO/Additional Patients/5.GSE46517/annotLookup.RData")
#reorder
indicesLookup <- match(rownames(gse72267), annotLookup$affy_hg_u133_plus_2)
#control
dftmp <- data.frame(rownames(gse72267), annotLookup[indicesLookup, c("affy_hg_u133_plus_2", "external_gene_name")])
head(dftmp, 20)
table(dftmp[,1] == dftmp[,2])
#dftmp <- dftmp[complete.cases(dftmp), ]

#cheat 
rownames(gse72267) <- paste(annotLookup[indicesLookup, "external_gene_name"], c(1:length(indicesLookup)), sep="_")
#gsub("_[0-9]*$", "", rownames(gse72267))

####Remove NA
rm <- grep("^NA_", rownames(gse72267))
gse72267 <- gse72267[-rm,]

#
gse72267 <- data.frame(t(gse72267))
gse72267$group <- rownames(gse72267)
blood_pd <- grep("^Blood_PD_", gse72267$group)
blood_ht <- grep("^Blood_healthy_", gse72267$group)

blood_pd <- cbind.data.frame(name='blood_pd',index=blood_pd, class= 1)
blood_ht <- cbind.data.frame(name='blood_ht',index=blood_ht, class= 0)
df <- rbind.data.frame(blood_ht, blood_pd)
newdata <- df[order(df$index),]

gse72267 <- cbind.data.frame(gse72267, newdata)
gse72267$group <- NULL
gse72267$name <- NULL
gse72267$index <- NULL

gse72267_ht <- subset(gse72267, class==0)
gse72267_pd <-  subset (gse72267, class==1)
