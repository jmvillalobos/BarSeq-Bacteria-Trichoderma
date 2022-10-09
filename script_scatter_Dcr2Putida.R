#analysis effect exudates of Dcr2 over psuedomonas putida
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_putida.txt", header=TRUE, row.names = 1, sep="\t")
colnames(heat)

heat[,"SumSpend"]= rowSums(heat[,c("M","N", "L")])/3
heat[,"control"]=rowSums(heat[,c("D", "F", "G", "H")])/4



genesFitdownDcr2=rownames(heat)[heat$SumSpend < -1 & heat$control > -1]
genesFitupDcr2=rownames(heat)[heat$SumSpend > 1 & heat$control < 1]
genesFitComple=rownames(heat)[heat$SumSpend < 1 & heat$control < -1.5 & heat$SumSpend  > -1]

length(genesFitdownDcr2)
length(genesFitupDcr2)
length(genesFitComple)

fitnenegWT=heat[genesFitdownDcr2,]
fitneposiWT=heat[genesFitupDcr2,]
fitneComple=heat[genesFitComple,]



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



#general scatter plot
plot( heat$SumSpend, heat$control,  col="darkgray", xlim = c(-5,2), ylim = c(-5,1), cex = 0.7,pch=19, 
      xlab = "T. atroviride Dcr2-Spent media", ylab= "Without Fungus", main = "Pseudomonas putida")
points(fitnenegWT$SumSpend, fitnenegWT$control, col="purple", pch=19,cex = 0.7  )
points(fitneposiWT$SumSpend, fitneposiWT$control, col="red", pch=19,cex = 0.7 )
points(fitneComple$SumSpend, fitneComple$control, col="blue", pch=19, cex = 0.7  )
#points for annotation
points(fitnenegNamesFA$SumSpend,fitnenegNamesFA$control, col="coral3", pch=19, cex = 1.2 )
points(fitnenegNamesABC$SumSpend,fitnenegNamesABC$control, col="green", pch=19, cex = 1.2 )
points(fitnenegNamesIron$SumSpend,fitnenegNamesIron$control, col="brown", pch=19, cex = 1.2 )

#text for ids
text(fitnenegWT$SumSpend, fitnenegWT$control, labels = row.names(fitnenegWT),col="black", cex=0.5, pos=3 )
text(fitneposiWT$SumSpend,fitneposiWT$control,labels = row.names(fitneposiWT),col="black", cex=0.5, pos=3 )
text(fitnenegNamesFA$SumSpend,fitnenegNamesFA$control,labels = row.names(fitnenegNamesFA),col="black", cex=0.5, pos=3 )
text(fitnenegNamesABC$SumSpend,fitnenegNamesABC$control, labels = row.names(fitnenegNamesABC),col="black", cex=0.5, pos=3 )
text(fitnenegNamesIron$SumSpend,fitnenegNamesIron$control, labels = row.names(fitnenegNamesIron),col="black", cex=0.5, pos=3 )
text(fitneComple$SumSpend,fitneComple$control, labels = row.names(fitneComple),col="black", cex=0.5, pos=3 )

legend(-3, -3, legend=c("Negative Fitness","Recovered","Fusaric acid resistance", "ABC", "Iron"),
       col=c("purple","blue","coral3", "green", "brown"), pch=19, cex=1.2 )

abline(h=0, col="black", lty=2, lwd=1)
abline(v=0, col="black", lty=2, lwd=1)
abline(-0.01,1)



Exporter as a PDF using 7X7

