###comparación pareada simiae

heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/28agost2022BarSeq/fitness_poly_fusputida.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"PolyPutida0.2"]= rowSums(heat["D1"])
heat[,"control"]=rowSums(heat[,c("A1","A2","A3")])/3

heat2=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_putida.txt", header=TRUE, row.names = 1, sep="\t")

heat2[,"SumSpend"]= rowSums(heat2[,c("A","B", "C")])/3
heat2[,"control"]=rowSums(heat2[,c("D", "F", "G", "H")])/4


genesFitdownWT=rownames(heat2)[heat2$SumSpend < -1 & heat2$control > -1]

PolyPutida2Neg=rownames(heat)[heat$PolyPutida0.2 < -2.5 & heat$control > -1]
PolyPutida2Up=rownames(heat)[heat$PolyPutida0.2 > 1 & heat$control < 1]
genesComple=rownames(heat)[heat$PolyPutida0.2 < 1 & heat$control < -1.5 & heat$PolyPutida0.2   > -1]

length(PolyPutida2Neg)
length(PolyPutida2Up)
length(genesComple)

fitneneg=heat[PolyPutida2Neg,]
fitneposi=heat[PolyPutida2Up,]
fitneComple=heat[genesComple,]

fitneMSessencials=heat[genesFitdownWT,]



#general scatter plot
plot(heat$PolyPutida0.2, heat$control,  col="darkgray", xlim = c(-6,1), ylim = c(-6,1), cex = 0.5,pch=19, 
     xlab =   "Polymyxin B 0.002mg/ml", ylab= "R2A media", main = "Pseudomonas putida KT2440")
points(fitneMSessencials$PolyPutida0.2,fitneMSessencials$control, col="green", pch=19, cex = 0.5 )


text(fitneneg$PolyPutida0.2, fitneneg$control,labels = row.names(fitneneg),col="black", cex=0.5, pos=3 )

fitneneg

legend(-1.5,-3, legend=c("Unaffected","SM Essentials"),
       col=c("darkgrey", "green"), pch=19, cex=1)
abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)


Exporter as a PDF using 7X7
