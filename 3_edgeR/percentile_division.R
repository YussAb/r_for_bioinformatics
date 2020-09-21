htseqcounts<- t(htseqcounts)
htseqcounts <- as.data.frame(htseqcounts)

Percentile_00  = quantile(htseqcounts$ENSG00000164327,0)
Percentile_25  = quantile(htseqcounts$ENSG00000164327, 0.25)
Percentile_75  = quantile(htseqcounts$ENSG00000164327, 0.75)
Percentile_100 = quantile(htseqcounts$ENSG00000164327)


RB = rbind(Percentile_00, Percentile_25, Percentile_75, Percentile_100)

dimnames(RB)[[2]] = "Value"

RB

htseqcounts$Group[htseqcounts$ENSG00000164327 >= Percentile_00 & htseqcounts$ENSG00000164327 <  Percentile_25]  = "high"
htseqcounts$Group[htseqcounts$ENSG00000164327 >= Percentile_25 & htseqcounts$ENSG00000164327 <  Percentile_75]  = "middle"
htseqcounts$Group[htseqcounts$ENSG00000164327 >= Percentile_75 & htseqcounts$ENSG00000164327 <= Percentile_100] = "low"

htseqcounts$Group = factor(htseqcounts$Group,
                    levels=c("high", "middle", "low"))
