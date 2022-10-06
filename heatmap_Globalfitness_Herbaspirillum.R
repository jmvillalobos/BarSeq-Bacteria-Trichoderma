########Heatmap to show fitness effect using the 3 strains of Trichoderma over Herbaspirillum

herbas_matrix=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_herbas.txt", header=TRUE, row.names = 1, sep="\t")

herbas_matrix[,"SumSpendWT"]= rowSums(herbas_matrix[,c("A","B", "C")])/3
herbas_matrix[,"SumSpendDcr2"]= rowSums(herbas_matrix[,c("M","N", "O")])/3
herbas_matrix[,"SumSpendTmk3"]= rowSums(herbas_matrix[,c("S","T", "U")])/3
herbas_matrix[,"control"]=rowSums(herbas_matrix[,c("D","E","F")])/3


genesFitdownWT=rownames(herbas_matrix)[herbas_matrix$SumSpendWT < -1 & herbas_matrix$control > -1]
signi=genesFitdownWT


herbas=herbas_matrix[signi,c("control","SumSpendWT","SumSpendTmk3")]



Herb=herbas[c(-23,-24,-25,-26,-27,-28,-30, -32),]

#Saving the matrix of fitness
write.table(Herb, file="/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_means_herb.txt",sep="\t")
# load package
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendextend)
#create heatmap using pheatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(Herb, 1, cal_z_score))
my_hclust_gene <- hclust(dist(data_subset_norm), method = "complete")

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 6)

clust=as.data.frame(my_gene_col)
l <- split(clust, clust$my_gene_col)
for(cl in names(l)) {
  write.table(l[[cl]], file=paste(cl, ".txt", sep=""), sep="\t")
  print(cl)
  print(length(l[[cl]]$my_gene_col))
}
###cortando arbol y guardando pdf

pheatmap(data_subset_norm, cluster_cols = FALSE, show_rownames = TRUE,
         cutree_rows = 3, color=hcl.colors(50, "BluYl"), cex=1) 

