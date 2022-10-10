###Comparisson R2A media and treatment with Polymyxin B 0.001mg/ml over P. putida
setwd("/Users/manuel/Documents/posdoc_berkeley/barseq_results/28agost2022BarSeq/")
###comparación pareada simiae

heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/28agost2022BarSeq/fitness_poly_fusputida.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"PolyPutida1"]= rowSums(heat[,c("C1","C2","C3","C4")])/4
heat[,"control"]=rowSums(heat[,c("A1","A2","A3")])/3


heat2=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_putida.txt", header=TRUE, row.names = 1, sep="\t")

heat2[,"SumSpend"]= rowSums(heat2[,c("A","B", "C")])/3
heat2[,"control"]=rowSums(heat2[,c("D", "F", "G", "H")])/4


genesFitdownWT=rownames(heat2)[heat2$SumSpend < -1 & heat2$control > -1]

PolyPutida1Neg=rownames(heat)[heat$PolyPutida1 < -1 & heat$control > -1]
PolyPutida1Up=rownames(heat)[heat$PolyPutida1 > 1 & heat$control < 1]
genesComple=rownames(heat)[heat$SumSpendWT < 1 & heat$control < -1.5 & heat$SumSpendWT  > -1]


length(PolyPutida1Neg)
length(PolyPutida1Up)
length(genesComple)

fitneneg=heat[PolyPutida1Neg,]
fitneposi=heat[PolyPutida1Up,]
fitneComple=heat[genesComple,]

fitneMSessencials=heat[genesFitdownWT,]

genes_interIron=c("PP_4755",
                  "PP_1224",
                  "PP_1898",
                  "PP_2378",
                  "PP_1082",
                  "PP_1083")

genes_interABC=c( 
  "PP_5327",
  "PP_5328",
  "PP_1779",
  "PP_1778"
)



genes_interFA=c("PP_3425",
                "PP_3426",
                "PP_3427",
                "PP_3428")



fitnenegNamesIron=heat[genes_interIron,]
fitnenegNamesFA=heat[genes_interFA,]
fitnenegNamesABC=heat[genes_interABC,]

fitneneg=heat[PolyPutida1Neg,]
fitneposi=heat[PolyPutida1Up,]
fitneComple=heat[genesComple,]
#general scatter plot
plot(heat$PolyPutida1, heat$control,  col="darkgray", xlim = c(-4,2), ylim = c(-4,2), cex = 0.5,pch=19, 
      xlab =   "Polymyxin B 0.001mg/ml", ylab= "R2A media", main = "Pseudomonas putida KT2440")
points(fitneMSessencials$PolyPutida1,fitneMSessencials$control, col="green", pch=19, cex = 0.5 )


text(fitneneg$PolyPutida1, fitneneg$control,labels = row.names(fitneneg),col="black", cex=0.5, pos=3 )


legend(-1.5,-5, legend=c("Unaffected","SM Essentials"),
       col=c("darkgrey", "green"), pch=19, cex=1)
abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)


Exporter as a PDF using 7X7
