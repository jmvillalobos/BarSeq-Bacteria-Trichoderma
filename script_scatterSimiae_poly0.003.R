###Comparisson R2A media and treatment with Polymyxin B 0.003mg/ml
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/28agost2022BarSeq/fitness_simiaePol_Fus.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"Polymyxn0.3"]= rowSums(heat[,c("C1","C2","C3","C4","C6")])/5
heat[,"control"]=rowSums(heat[,c("A1","A2")])/2


heatWT=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_simiae.txt", header=TRUE, row.names = 1, sep="\t")
heatWT[,"SumSpend"]= rowSums(heatWT[,c("A","B","C")])/3
heatWT[,"control"]=rowSums(heatWT[,c("D","E")])/2


genesFitdownWT=rownames(heatWT)[heatWT$SumSpend < -1 & heatWT$control > -1]


genesFitdownWTpoly3=rownames(heat)[heat$Polymyxn0.3 < -1 & heat$control > -1]
genesFitupWTpoly3=rownames(heat)[heat$Polymyxn0.3 > 1 & heat$control < 1]
genesFitComple=rownames(heat)[heat$Polymyxn0.3 < 1 & heat$control < -1.5& heat$Polymyxn0.3 > -1]


#genesFitnegCon=rownames(heat)[heat$SumSpendWT < 1 & heat$control < -1.5& heat$SumSpendWT  > -1]
essencialSM=heat[genesFitdownWT,]

fitnenegPoly3=heat[genesFitdownWTpoly3,]
fitneposiPoly3=heat[genesFitupWTpoly3,]
fitComplePoly3=heat[genesFitComple,]


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

#plot of the gene fitness in response to pure polymyxin B.

plot( heat$Polymyxn0.3, heat$control,  col="darkgray", xlim = c(-8,3), ylim = c(-8,2), cex = 0.5,pch=19,
      xlab = "Polymyxin B 0.003mg/ml", ylab= "R2A media", main = "Pseudomonas simiae WCS417")
points(essencialSM$Polymyxn0.3,essencialSM$control, col="green", pch=19, cex = 0.5 )
points(fitnenegNamesArn$Polymyxn0.3,fitnenegNamesArn$control, col="brown", pch=19, cex = 1 )
text(fitnenegNamesArn$Polymyxn0.3,fitnenegNamesArn$control,labels = row.names(fitnenegNamesArn),col="black", cex=0.5, pos=3 )
legend(-1.5,-5, legend=c("SM Essentials", "Arn-operon"),
       col=c("green", "brown"), pch=19, cex=1)
abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)
