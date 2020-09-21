#df = data frame con gene expression
#x = gene of interest

Percentile_00  = quantile(df$x,0)
Percentile_25  = quantile(df$x, 0.25)
Percentile_75  = quantile(df$x, 0.75)
Percentile_100 = quantile(df$x,1)
#RB = rbind(Percentile_00, Percentile_25, Percentile_75, Percentile_100)
#dimnames(RB)[[2]] = "Value"
#RB

df$Group[df$x >= Percentile_00 & df$x <  Percentile_25]  = "low"
df$Group[df$x >= Percentile_25 & df$x <  Percentile_75]  = "middle"
df$Group[df$x >= Percentile_75 & df$x <= Percentile_100] = "high"
dfkm <-df[!(df$Group == "middle"), ]


#SURVIVAL ANALYSIS 1
#Kaplan-Meier non-parametric analysis

kmsurvival1 <- survfit(Surv(time, event)~ Group, data=PT)
ggsurvplot(kmsurvival1, linetype = "strata", 
           conf.int = TRUE, pval = TRUE,
           palette = "Dark2")

kmsurvival2 <- survfit(Surv(time, event)~ Group, data=PT)
ggsurvplot(kmsurvival2, linetype = "strata", 
           conf.int = TRUE, pval = TRUE,
           palette = "Dark2")

kmsurvival3 <- survfit(Surv(time, event)~ Group, data=MT)
ggsurvplot(kmsurvival3, linetype = "strata", 
           conf.int = TRUE, pval = TRUE,
           palette = "Dark2")

kmsurvival4 <- survfit(Surv(time, event)~ Group, data=MTkm)
ggsurvplot(kmsurvival4, linetype = "strata", 
           conf.int = TRUE, pval = TRUE,
           palette = "Dark2")

