#TGCABiolinks
#sotto di forma conte
#HTseq counts
#no fpkm ecc.
##https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("", version = "3.8")

library("TCGAbiolinks")
library("DT")

DataDirectory <- paste0("../GDC/",gsub("-","_","TCGA-SKCM"))
FileNameData <- paste0(DataDirectory, "_","HTSeq_Counts",".rda")

query <- GDCquery(project = "TCGA-SKCM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols = c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, #TP primary solid tumor
                                  typesample = "TP")

dataSmTM <- TCGAquery_SampleTypes(barcode = samplesDown, 
                                  typesample = "TM")

queryDown <- GDCquery(project = "TCGA-SKCM", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTM))

GDCdownload(query = queryDown)

dataPrep1 <- GDCprepare(query = queryDown, 
                        save = TRUE, 
                        save.filename = "TCGA_SKCM_HTSeq_Countds.rda")

#############################################################################
library(SummarizedExperiment)

datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

pheno <- as.data.frame(colData(dataPrep1))

datatable(assay(data)[1:100,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

htseqcounts <- assay(dataPrep1)

x <-htseqcounts["ENSG00000108064",] #tfam
y<- htseqcounts["ENSG00000164327",] #rictor

