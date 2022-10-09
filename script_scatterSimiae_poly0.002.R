###comparisson 0.002mg/ml polymyxin B.

heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/28agost2022BarSeq/fitness_simiaePol_Fus.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"Polymyxn0.2"]= rowSums(heat[,c("B1","B2","B3","B4","B5","B6")])/6
heat[,"control"]=rowSums(heat[,c("A1","A2")])/2



genesFitdownWTpoly2=rownames(heat)[heat$Polymyxn0.2 < -1 & heat$control > -1]
genesFitupWTpoly2=rownames(heat)[heat$Polymyxn0.2 > 1 & heat$control < 1]
genesFitComple=rownames(heat)[heat$Polymyxn0.2 < 1 & heat$control < -1.5& heat$Polymyxn0.3 > -1]

heatWT=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_simiae.txt", header=TRUE, row.names = 1, sep="\t")
heatWT[,"SumSpend"]= rowSums(heatWT[,c("A","B","C")])/3
heatWT[,"control"]=rowSums(heatWT[,c("D","E")])/2


genesFitdownWT=rownames(heatWT)[heatWT$SumSpend < -1 & heatWT$control > -1]

essencialSM=heat[genesFitdownWT,]
fitnenegPoly2=heat[genesFitdownWTpoly2,]
fitneposiPoly2=heat[genesFitupWTpoly2,]
fitComplePoly2=heat[genesFitComple,]


genes_interArn=c("PS417_13800",
                 "PS417_13805",
                 "PS417_13795",
                 "PS417_13790")

fitnenegNamesArn=heat[genes_interArn,]


#plot of the gene fitness in response to pure polymyxin B.

plot( heat$Polymyxn0.2, heat$control,  col="darkgray", xlim = c(-5,3), ylim = c(-5,2), cex = 0.5,pch=19,
      xlab = "Polymyxin B 0.002mg/ml", ylab= "R2A media", main = "Pseudomonas simiae WCS417")
points(essencialSM$Polymyxn0.2,essencialSM$control, col="green", pch=19, cex = 0.5 )
points(fitnenegNamesArn$Polymyxn0.2,fitnenegNamesArn$control, col="brown", pch=19, cex = 1 )
text(essencialSM$Polymyxn0.2,essencialSM$control,labels = row.names(essencialSM),col="black", cex=0.5, pos=3 )
text(fitnenegNamesArn$Polymyxn0.2,fitnenegNamesArn$control,labels = row.names(fitnenegNamesArn),col="black", cex=0.5, pos=3 )
legend(-1.5,-3, legend=c("SM Essentials", "Arn-operon"),
       col=c("green", "brown"), pch=19, cex=1)
abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)


