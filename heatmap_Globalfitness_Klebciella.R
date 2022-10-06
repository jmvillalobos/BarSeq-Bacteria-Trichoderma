####Heatmap to show fitness effect using the 3 strains of Trichoderma.
kleb_matrix=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_klebsiella.txt", header=TRUE, row.names = 1, sep="\t")

kleb_matrix[,"SumSpendWT"]= rowSums(kleb_matrix[,c("A","B","C")])/3
kleb_matrix[,"SumSpendDcr2"]= rowSums(kleb_matrix[,c("K","L","M")])/3
kleb_matrix[,"SumSpendTmk3"]= rowSums(kleb_matrix[,c("Q","R","S")])/3
kleb_matrix[,"control"]=rowSums(kleb_matrix[,c("D","F","G")])/3

genesFitdownWT=rownames(kleb_matrix)[kleb_matrix$SumSpendWT < -1.5 & kleb_matrix$control > -1]
signi=genesFitdownWT


kleb=kleb_matrix[signi,c("control","SumSpendWT","SumSpendTmk3")]

#saving the table for the heatmap with genes with significants fitnees in the assay.
#write.table(kleb, file="/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_means_kleb.txt",sep="\t")


##for  the heatmap
# load package
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendextend)
#create heatmap using pheatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(kleb, 1, cal_z_score))
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
         cutree_rows = 3, color=hcl.colors(50, "BluYl"), cex=1) 
