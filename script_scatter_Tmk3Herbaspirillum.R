###tmk3
#analysis psuedomonas
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_herbas.txt", header=TRUE, row.names = 1, sep="\t")
colnames(heat)
heat[,"SumSpend"]= rowSums(heat[,c("S","T", "U")])/3
heat[,"control"]=rowSums(heat[,c("D","E","F")])/3


genesFitdown=rownames(heat)[heat$SumSpend < -1 & heat$control > -1]
genesFitup=rownames(heat)[heat$SumSpend > 1 & heat$control < 1]
=rownames(heat)[heat$SumSpend < 1 & heat$control < -2& heat$SumSpend  > -2]


length(genesFitdown)
length(genesFitup)
length(genesFitComple)

fitneneg=heat[genesFitdown,]
fitneposi=heat[genesFitup,]
fitneComple=heat[genesFitComple,]



comple=c( "HSERO_RS05910",
          "HSERO_RS07625",
          "HSERO_RS15630",
          "HSERO_RS16070",
          "HSERO_RS19995",
          "HSERO_RS21035")



genes_interIron=c("HSERO_RS02660",
                  "HSERO_RS03450",
                  "HSERO_RS00165",
                  "HSERO_RS03440",
                  "HSERO_RS13650",
                  "HSERO_RS03445")



genes_interABC=c( 
  "HSERO_RS10850"
)



genes_interArn=c("HSERO_RS21375")



fitnenegNamesIron=heat[genes_interIron,]
fitnenegNamesABC=heat[genes_interABC,]
fitnenegNamesArn=heat[genes_interArn,]



#general scatter plot
plot( heat$SumSpend, heat$control,  col="darkgray", xlim = c(-4,1), ylim = c(-4,1), cex = 0.7,pch=19, 
      xlab = "T. atroviride Tmk3-Spent media", ylab= "Without Fungus", main = "Herbaspirillum seropedicae SmR1")
points(fitneneg$SumSpend, fitneneg$control, col="purple", pch=19,cex = 0.7  )
points(fitneposi$SumSpend, fitneposi$control, col="red", pch=19,cex = 0.7 )
points(fitneComple$SumSpend, fitneComple$control, col="blue", pch=19, cex = 0.7  )
#points for annotation
points(fitnenegNamesIron$SumSpend,fitnenegNamesIron$control, col="coral3", pch=19, cex = 1.2 )
points(fitnenegNamesABC$SumSpend,fitnenegNamesABC$control, col="green", pch=19, cex = 1.2 )
points(fitnenegNamesArn$SumSpend,fitnenegNamesArn$control, col="brown", pch=19, cex = 1.2 )

#text for ids
text(fitneneg$SumSpend, fitneneg$control, labels = row.names(fitneneg),col="black", cex=0.5, pos=3 )

text(fitnenegNamesIron$SumSpend,fitnenegNamesIron$control,labels = row.names(fitnenegNamesIron),col="black", cex=0.5, pos=3 )
text(fitnenegNamesABC$SumSpend,fitnenegNamesABC$control, labels = row.names(fitnenegNamesABC),col="black", cex=0.5, pos=3 )
text(fitnenegNamesArn$SumSpend,fitnenegNamesArn$control, labels = row.names(fitnenegNamesArn),col="black", cex=0.5, pos=3 )
text(fitneComple$SumSpend,fitneComple$control, labels = row.names(fitneComple),col="black", cex=0.5, pos=3 )

legend(-1.5, -3, legend=c("Negative Fitness","Recovered","Iron", "ABC", "Arn-operon"),
       col=c("purple","blue","coral3", "green", "brown"), pch=19, cex=1.2 )

abline(h=0, col="black", lty=2, lwd=1)
abline(v=0, col="black", lty=2, lwd=1)
abline(-0.01,1)



Exporter as a PDF using 7X7