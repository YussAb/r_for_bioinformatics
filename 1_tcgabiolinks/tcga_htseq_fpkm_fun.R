#https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/250

tcga_htseq <- function(x){
library('TCGAbiolinks') 
DataDirectory <- paste0("../GDC/",gsub("-","_", paste("TCGA-", x, sep= ""))) 
FileNameData <- paste0(DataDirectory, "_","HTSeq_FPKM",".rda")

query <- GDCquery(project = paste("TCGA-", x, sep= ""),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")

samplesDown <- getResults(query,cols = c("cases"))

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown, #TP primary solid tumor
                                  typesample = "TP")

dataSmTM <- TCGAquery_SampleTypes(barcode = samplesDown, #TM metastases
                                  typesample = "TM")

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown, #NT normal tissue
                                  typesample = "NT")

queryDown <- GDCquery(project = paste("TCGA-", x, sep= ""), 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - FPKM", 
                      barcode = c(dataSmTP, dataSmTM,dataSmNT))

GDCdownload(query = queryDown)

dataPrep1 <- GDCprepare(query = queryDown, 
                        save = TRUE, 
                        save.filename = paste ("TCGA_", x, "_HTSeq_FPKM.rda", sep=""))
return(dataPrep1)
}
data <- tcga_htseq('UVM')

#############################################################################
library(SummarizedExperiment)

pheno <- as.data.frame(colData(data))

htseqfpkm <- assay(data)

x<-htseqfpkm["ENSG00000108064",] #tfam
y<- htseqfpkm["ENSG00000164327",] #rictor



