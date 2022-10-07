#analysis dcr2 over psuedomonas
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_simiae.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"SumSpend"]= rowSums(heat["J"])
heat[,"control"]=rowSums(heat[,c("D","E")])/2

genesFitdown=rownames(heat)[heat$SumSpend < -1.5 & heat$control > -1]
genesFitup=rownames(heat)[heat$SumSpend > 1.5 & heat$control < 1]
genesFitComple=rownames(heat)[heat$SumSpend < 1 & heat$control < -1.5& heat$SumSpend  > -1]

length(genesFitdown)
length(genesFitup)
length(genesFitComple)

fitneneg=heat[genesFitdown,]
fitneposi=heat[genesFitup,]
fitneComple=heat[genesFitComple,]


genes_interFA=c("PS417_04740",
                "PS417_04745",
                "PS417_04750",
                "PS417_04755")

genes_interABC=c("PS417_10720",
                 "PS417_24530",
                 "PS417_04380")


genes_interArn=c("PS417_13790",
                 "PS417_13795",
                 "PS417_13800",
                 "PS417_13805",
                 "PS417_13810",
                 "PS417_13815",
                 "PS417_13820")



fitnenegNamesFA=heat[genes_interFA,]
fitnenegNamesABC=heat[genes_interABC,]
fitnenegNamesArn=heat[genes_interArn,]


#general scatter plot
plot( heat$SumSpend, heat$control,  col="darkgray", xlim = c(-4,2), ylim = c(-4,1), cex = 0.7,pch=19, 
      xlab = "T. atroviride Dcr2-Spent media", ylab= "Without Fungus", main = "Pseudomonas simiae WCS417")
points(fitneneg$SumSpend, fitneneg$control, col="purple", pch=19,cex = 0.7  )
points(fitneposi$SumSpend, fitneposi$control, col="red", pch=19,cex = 0.7 )
points(fitneComple$SumSpend, fitneComple$control, col="blue", pch=19, cex = 0.7  )
#points for annotation
points(fitnenegNamesFA$SumSpend,fitnenegNamesFA$control, col="coral3", pch=19, cex = 1.2 )
points(fitnenegNamesABC$SumSpend,fitnenegNamesABC$control, col="green", pch=19, cex = 1.2 )
points(fitnenegNamesArn$SumSpend,fitnenegNamesArn$control, col="brown", pch=19, cex = 1.2 )

#text for ids
text(fitneneg$SumSpend, fitneneg$control, labels = row.names(fitneneg),col="black", cex=0.5, pos=3 )

text(fitnenegNamesFA$SumSpend,fitnenegNamesFA$control,labels = row.names(fitnenegNamesFA),col="black", cex=0.5, pos=3 )
text(fitnenegNamesABC$SumSpend,fitnenegNamesABC$control, labels = row.names(fitnenegNamesABC),col="black", cex=0.5, pos=3 )
text(fitnenegNamesArn$SumSpend,fitnenegNamesArn$control, labels = row.names(fitnenegNamesArn),col="black", cex=0.5, pos=3 )
text(fitneComple$SumSpend,fitneComple$control, labels = row.names(fitneComple),col="black", cex=0.5, pos=3 )

legend(0.5, -3, legend=c("Negative Fitness","Recovered","Fusaric acid resistance", "ABC", "Arn-operon"),
       col=c("purple","blue","coral3", "green", "brown"), pch=19, cex=1.2 )

abline(h=0, col="black", lty=2, lwd=1)
abline(v=0, col="black", lty=2, lwd=1)
abline(-0.01,1)



Exporter as a PDF using 7X7

