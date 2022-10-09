###Comparisson R2A media and treatment with Fusaric Acid 120 ug/ml
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/28agost2022BarSeq/fitness_simiaePol_Fus.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"Fus120"]= rowSums(heat[,c("D1","D2","D3","D4","D6")])/6
heat[,"control"]=heat$E1


heatWT=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_simiae.txt", header=TRUE, row.names = 1, sep="\t")
heatWT[,"SumSpend"]= rowSums(heatWT[,c("A","B","C")])/3
heatWT[,"control"]=rowSums(heatWT[,c("D","E")])/2

genesFitdownWT=rownames(heatWT)[heatWT$SumSpend < -1 & heatWT$control > -1]

genesFitdownWTFus=rownames(heat)[heat$Fus120 < -1 & heat$control > -1]
genesFitupWTFus=rownames(heat)[heat$Fus120 > 1 & heat$control < 1]
genesFitComple=rownames(heat)[heat$Fus120 < 1 & heat$control < -1.5& heat$Fus120 > -1]

length(genesFitdownWTFus)
length(genesFitupWTFus)
length(genesFitComple)
#genesFitnegCon=rownames(heat)[heat$SumSpendWT < 1 & heat$control < -1.5& heat$SumSpendWT  > -1]
essencialSM=heat[genesFitdownWT,]

fitnenegFus=heat[genesFitdownWTFus,]
fitneposiFus=heat[genesFitupWTFus,]
fitCompleFus=heat[genesFitComple,]

genes_interFA=c("PS417_04740",
                "PS417_04745",
                "PS417_13190")


genes_interABC=c( 
  "PS417_10720",
  "PS417_24530",
  "PS417_04380")


fitnenegNamesFA=heat[genes_interFA,]
fitnenegNamesABC=heat[genes_interABC,]


#plot of the gene fitness in response to pure polymyxin B.

plot( heat$Fus120, heat$control,  col="darkgray", xlim = c(-5,3), ylim = c(-5,2), cex = 0.5,pch=19,
      xlab = "Fusaric Acid 120 ug/ml", ylab= "R2A media", main = "Pseudomonas simiae WCS417")
points(essencialSM$Fus120,essencialSM$control, col="green", pch=19, cex = 0.5 )
points(fitnenegNamesFA$Fus120,fitnenegNamesFA$control, col="brown", pch=19, cex = 1 )
points(fitnenegNamesABC$Fus120,fitnenegNamesABC$control, col="darkgreen", pch=19, cex = 1 )

text(fitnenegFus$Fus120,fitnenegFus$control,labels = row.names(fitnenegFus),col="black", cex=0.5, pos=3 )
text(fitnenegNamesABC$Fus120,fitnenegNamesABC$control,labels = row.names(fitnenegNamesABC),col="black", cex=0.5, pos=3 )
text(fitnenegNamesFA$Fus120,fitnenegNamesFA$control,labels = row.names(fitnenegNamesFA),col="black", cex=0.5, pos=3 )
legend(-2,-3, legend=c("SM Essentials", "Fusaric acid resistance", "ABC transporters"),
       col=c("green", "brown", "darkgreen"), pch=19, cex=1)
abline(h=-1, col="black", lty=2, lwd=1)
abline(v=-1.5, col="black", lty=2, lwd=1)
abline(-0.01,1)
