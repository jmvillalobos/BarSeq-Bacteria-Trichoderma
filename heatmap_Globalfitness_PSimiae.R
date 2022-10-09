###heatmaps
simiae_matrix=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_simiae.txt", header=TRUE, row.names = 1, sep="\t")
simiae_matrix[,"SumSpend"]= rowSums(simiae_matrix[,c("A","B","C")])/3
simiae_matrix[,"control"]=rowSums(simiae_matrix[,c("D","E")])/2
simiae_matrix[,"SumSpendDcr2"]= rowSums(simiae_matrix["J"])
simiae_matrix[,"SumSpendtmk3"]= rowSums(simiae_matrix[c("L", "O")])/2

genesFitdownWT=rownames(simiae_matrix)[simiae_matrix$SumSpend < -1 & simiae_matrix$control > -1]
signi=genesFitdownWT


Simiae=simiae_matrix[signi,c("control","SumSpend","SumSpendtmk3")]

#saving the table for the heatmap with genes with significants fitnees in the assay.
#write.table(Simiae, file="/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_means_simiae.txt",sep="\t")


# load package
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendextend)
#create heatmap using pheatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(Simiae, 1, cal_z_score))
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
#pdf("heatmap_cpm_cutree_6.pdf") 
pheatmap(data_subset_norm, cluster_cols = FALSE, show_rownames = TRUE,
         cutree_rows = 3, color=hcl.colors(50, "BluYl"), cex=0.8
) 

