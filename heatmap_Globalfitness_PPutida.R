####heatmap for putida and mutants

putida_matrix=read.table("/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_putida.txt", header=TRUE, row.names = 1, sep="\t")

putida_matrix[,"SumSpendWT"]= rowSums(putida_matrix[,c("A","B", "C")])/3
putida_matrix[,"SumSpendDcr2"]= rowSums(putida_matrix[,c("M","N", "L")])/3
putida_matrix[,"SumSpendTmk3"]= rowSums(putida_matrix[,c("R","S", "T")])/3
putida_matrix[,"control"]=rowSums(putida_matrix[,c("D", "F", "G", "H")])/4

genesFitdownWT=rownames(putida_matrix)[putida_matrix$SumSpendWT< -1 & putida_matrix$control > -1]

signi=genesFitdownWT

putida=putida_matrix[signi,c("control","SumSpendWT","SumSpendTmk3")]

putida=putida[c(-15,-16,-17,-18,-33,-34,-38,-40,-41,-42,-44,-45,-46),]

####saving table fo the fitness in putida
#write.table(putida, file="/Users/manuel/Documents/posdoc_berkeley/barseq_results/tabla_means_putida.txt",sep="\t")

# load package
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendextend)
#create heatmap using pheatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
data_subset_norm <- t(apply(putida, 1, cal_z_score))
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
         cutree_rows = 3, color=hcl.colors(50, "BluYl"), cex=1
) 
