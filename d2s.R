tree<-upgma(D2S_matrix, method="average")
plot(tree, main="UPGMA")

object_name <- as.dist(D2S_matrix)
clust.res<-hclust(object_name, method = "single")
clust.res


plot(clust.res, ylab = "D2S dissimilarity")
