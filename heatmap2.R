library(d3heatmap)
Location<- read.csv("C:/Users/Rachael Karns/Desktop/top_taxa.csv", header = T, row.names = 1)
Location<-Location[,order(colnames(Location))]

Location<-data.matrix(Location, rownames.force = NA)
Location<-t(Location)
d3heatmap(Location, Colv = NA,Rowv = NA, col = c("gray96", "turquoise4"), yaxis_width=290, cexRow = 0.6,cexCol = 1)

