#######################################################
#analysis psuedomonas putida KT2440
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_putida.txt", header=TRUE, row.names = 1, sep="\t")
colnames(heat)


heat[,"SumSpend"]= rowSums(heat[,c("A","B", "C")])/3
heat[,"control"]=rowSums(heat[,c("D", "F", "G", "H")])/4
genesFitWTdown=rownames(heat)[heat$SumSpend < -1& heat$control > -1]
genesFitWTup=rownames(heat)[heat$SumSpend > 1 & heat$control < 1]

###comparison with  Fusaric acid in differents concentrations
heat2=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/28agost2022BarSeq/fitness_poly_fusputida.txt", header=TRUE, row.names = 1, sep="\t")
heat2[,"Fusaric480"]= rowSums(heat2[c("G1","G2","G3","G4","G5")])/5
heat2[,"Fusaric240"]= rowSums(heat2[c("F1","F2","F3","F4","F5")])/5
heat2[,"Fusaric120"]= rowSums(heat2[c("E1","E2","E3","E4","E5")])/5
heat2[,"control"]=rowSums(heat2[,c("H1","H2","H3","H4","H5")])/5


fitnesupernegFus480=rownames(heat2)[heat2$Fusaric480 < -1 & heat2$control > -1]
fitnesupernegFus240=rownames(heat2)[heat2$Fusaric240 < -1 & heat2$control > -1]
fitnesupernegFus120=rownames(heat2)[heat2$Fusaric120 < -1 & heat2$control > -1]
negFus480=heat2[fitnesupernegFus480,]
negFus240=heat2[fitnesupernegFus240,]
negFus120=heat2[fitnesupernegFus120,]
fitnesnegWT=heat2[genesFitWTdown,]

###plot for P. putida responding to 480 ug/ml
plot( heat2$Fusaric480, heat2$control,  col="darkgray", xlim = c(-6,1), ylim = c(-6,1), cex = 0.5,pch=19,
      xlab =   "Fusaric acid 480 ug/ml", ylab= "R2A media", main = "Pseudomonas putida KT2440")
points(negFus480$Fusaric480, negFus480$control, col="purple", pch=19,cex = 0.5  )
points(fitnesnegWT$Fusaric480, fitnesnegWT$control, col="green", pch=19,cex = 0.8  )

text(negFus480$Fusaric480, negFus480$control, labels = row.names(negFus480),col="black", cex=0.5, pos=3)

abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)

legend(-2,-4, legend=c("Unaffected","Affected by FA","SM Essentials"),
       col=c("darkgray","purple", "green"), pch=19, cex=1)



###plot for P. putida responding to 240 ug/ml
plot( heat2$Fusaric240, heat2$control,  col="darkgray", xlim = c(-6,1), ylim = c(-6,1), cex = 0.5,pch=19,
      xlab =   "Fusaric acid 240 ug/ml", ylab= "R2A media", main = "Pseudomonas putida KT2440")
points(negFus240$Fusaric240, negFus240$control, col="purple", pch=19,cex = 0.5  )
points(fitnesnegWT$Fusaric240, fitnesnegWT$control, col="green", pch=19,cex = 0.8  )

text(negFus240$Fusaric240, negFus240$control, labels = row.names(negFus240),col="black", cex=0.5, pos=3)

abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)

legend(-2,-4, legend=c("Unaffected","Affected by FA","SM Essentials"),
       col=c("darkgray","purple", "green"), pch=19, cex=1)


###plot for P. putida responding to 120 ug/ml
plot( heat2$Fusaric120, heat2$control,  col="darkgray", xlim = c(-6,1), ylim = c(-6,1), cex = 0.5,pch=19,
      xlab =   "Fusaric acid 240 ug/ml", ylab= "R2A media", main = "Pseudomonas putida KT2440")
points(negFus120$Fusaric120, negFus120$control, col="purple", pch=19,cex = 0.5  )
points(fitnesnegWT$Fusaric120, fitnesnegWT$control, col="green", pch=19,cex = 0.8  )

text(negFus120$Fusaric120, negFus120$control, labels = row.names(negFus120),col="black", cex=0.5, pos=3)

abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)

legend(-2,-4, legend=c("Unaffected","Affected by FA","SM Essentials"),
       col=c("darkgray","purple", "green"), pch=19, cex=1)




