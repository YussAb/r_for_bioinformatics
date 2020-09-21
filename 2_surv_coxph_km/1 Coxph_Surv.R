#NB####pheno is from tcgabiolinks

###Time
#1
SKCMpheno <- read.delim("TCGA-SKCM.GDC_phenotype.tsv")
sample <- select(pheno,sample)
daysdeath <- select(pheno,days_to_death) #CX
dayslastfollow <- select(pheno, days_to_last_follow_up) #CY
time <- cbind.data.frame (daysdeath, dayslastfollow)
time$days_to_death[is.na(time$days_to_death)] <- time$days_to_last_follow_up[is.na(time$days_to_death)]
daysdeath <- time[[1]]


###Event
vitalstatus <- select(pheno, vital_status) #DI

vs <- function(x){
if (x == 'alive'){
  x <- 0
} else {
  x <-  1
}
return(x)}

vitalstatus$vs <- sapply(vitalstatus$vital_status, vs)

SDV <- data.frame(sample, vitalstatus$vs, daysdeath)
colnames(SDV) <- c("sample", "event", "time")

#Cox proportional hazard model - coefficients and hazard rates
#define variables
#time<- PT$time
#event<- PT$event
#Rictor <- PT$RICTOR
#mitoscore <- PT$MitoScore
library(survival)

coxph <- coxph(Surv(SDV$time, SDV$event) ~ htseqfpkm[1,])
coxph
