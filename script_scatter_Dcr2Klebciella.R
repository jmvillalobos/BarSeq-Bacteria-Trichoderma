###Scatter plot over the result with exudates of dcr2 mutant in Trichoderma

#analysis klebsiella
heat=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_klebsiella.txt", header=TRUE, row.names = 1, sep="\t")
heat[,"SumSpend"]= rowSums(heat[,c("K","L","M")])/3
heat[,"control"]=rowSums(heat[,c("D","F","G")])/3


genesFitdown=rownames(heat)[heat$SumSpend < -1.5 & heat$control > -1]
genesFitup=rownames(heat)[heat$SumSpend > 1.5 & heat$control < 1]
genesFitComple=rownames(heat)[heat$SumSpend < 1 & heat$control < -2& heat$SumSpend  > -2]


fitneneg=heat[genesFitdown,]
fitneposi=heat[genesFitup,]
fitneComple=heat[genesFitComple,]


###interesting genes
transcrip=c("BWI76_RS00760",
            "BWI76_RS05015",
            "BWI76_RS07465",
            "BWI76_RS09825",
            "BWI76_RS11275",
            "BWI76_RS20360",
            "BWI76_RS06765",
            "BWI76_RS16250",
            "BWI76_RS21310",
            "BWI76_RS07465")

genes_interIron=c("BWI76_RS17190",
                  "BWI76_RS24275",
                  "BWI76_RS11135",
                  "BWI76_RS24280",
                  "BWI76_RS05135",
                  "BWI76_RS08535",
                  "BWI76_RS08525",
                  "BWI76_RS26620"
)

genes_interABC=c( 
  "BWI76_RS09640",
  "BWI76_RS26395",
  "BWI76_RS26390",
  "BWI76_RS00760")

genes_interArn=c("BWI76_RS26460",
                 "BWI76_RS26455",
                 "BWI76_RS11275")

fitnenegNamesTra=heat[transcrip,]
fitnenegNamesIron=heat[genes_interIron,]
fitnenegNamesABC=heat[genes_interABC,]
fitnenegNamesArn=heat[genes_interArn,]


#general scatter plot
plot( heat$SumSpend, heat$control,  col="darkgray", xlim = c(-8,1), ylim = c(-8,1), cex = 0.7,pch=19, 
      xlab = "T. atroviride Dcr2-Spent media", ylab= "Without Fungus", main = "Klebsiella michiganensis M5aI")
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

legend(-3, -6, legend=c("Negative Fitness","Recovered","Iron", "ABC", "Arn-operon"),
       col=c("purple","blue","coral3", "green", "brown"), pch=19, cex=1.2 )

abline(h=0, col="black", lty=2, lwd=1)
abline(v=0, col="black", lty=2, lwd=1)
abline(-0.01,1)



Exporter as a PDF using 7X7
