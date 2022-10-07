#analysis psuedomonas
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_simiae.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"SumSpend"]= rowSums(heat2[,c("A","B","C")])/3
heat[,"control"]=rowSums(heat2[,c("D","E")])/2


genesFitdownWT=rownames(heat)[heat$SumSpend < -1.5 & heat$control > -1]
genesFitupWT=rownames(heat)[heat$SumSpend > 1.5 & heat$control < 1]
genesFitComple=rownames(heat)[heat$SumSpend < 1 & heat$control < -2& heat$SumSpend  > -2]


fitnenegWT=heat[genesFitdownWT,]
fitneposiWT=heat[genesFitupWT,]
fitneComple=heat[genesFitComple,]

genes_interFA=c("PS417_04740",
                "PS417_04745",
                "PS417_04750",
                "PS417_04755",
                "PS417_00710",
                "PS417_00715", 
                "PS417_00720",
                "PS417_00725")


genes_interABC=c("PS417_10720",
                 "PS417_24530",
                 "PS417_04380")


genes_interArn=c("PS417_13790",
                 "PS417_13800",
                 "PS417_13805",
                 "PS417_13795",
                 "PS417_13790")



fitnenegNamesFA=heat[genes_interFA,]
fitnenegNamesABC=heat[genes_interABC,]
fitnenegNamesArn=heat[genes_interArn,]


#general scatter plot
plot( heat$SumSpend, heat$control,  col="darkgray", xlim = c(-4,2), ylim = c(-4,1), cex = 0.7,pch=19, 
      xlab = "T. atroviride WT-Spent media", ylab= "Without Fungus", main = "Pseudomonas simiae WCS417")
points(fitnenegWT$SumSpend, fitnenegWT$control, col="purple", pch=19,cex = 0.7  )
points(fitneposiWT$SumSpend, fitneposiWT$control, col="red", pch=19,cex = 0.7 )
points(fitneComple$SumSpend, fitneComple$control, col="blue", pch=19, cex = 0.7  )
#points for annotation
points(fitnenegNamesFA$SumSpend,fitnenegNamesFA$control, col="coral3", pch=19, cex = 1.2 )
points(fitnenegNamesABC$SumSpend,fitnenegNamesABC$control, col="green", pch=19, cex = 1.2 )
points(fitnenegNamesArn$SumSpend,fitnenegNamesArn$control, col="brown", pch=19, cex = 1.2 )

#text for ids
text(fitnenegWT$SumSpend, fitnenegWT$control, labels = row.names(fitnenegWT),col="black", cex=0.5, pos=3 )

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

