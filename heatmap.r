#Heatmap generator
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: December 3, 2019 by Rachael Storo

#Purpose: Create a script for visualizing relative abundance tables as a heat map


library(d3heatmap)
library(htmlwidgets)
top<- read.csv("/work/sbjlab/rck/Sed_Assemblies/AT26-13-87/87_order_profile.csv", header = T, row.names = 1)
map1<-d3heatmap(top, Colv = NA,Rowv = NA, col = c("grey18", "forestgreen"), cexRow = 0.6,cexCol = 1)
saveWidget(map1, "87.html")

mid<- read.csv("/work/sbjlab/rck/Sed_Assemblies/AT26-13-89/89_order_profile.csv", header = T, row.names = 1)
map2<-d3heatmap(mid, Colv = NA,Rowv = NA, col = c("grey18", "forestgreen"), cexRow = 0.6,cexCol = 1)
saveWidget(map2, "89.html")

low<- read.csv("/work/sbjlab/rck/Sed_Assemblies/AT26-13-91/91_order_profile.csv", header = T, row.names = 1)
map3<-d3heatmap(low, Colv = NA,Rowv = NA, col = c("grey18", "forestgreen"), cexRow = 0.6,cexCol = 1)
saveWidget(map3, "91.html")

all<- read.csv("/work/sbjlab/rck/Sed_Assemblies/combined_87-89-91_Metaphlan_profile.csv", header = T, row.names = 1)
map4<-d3heatmap(all, Colv = NA,Rowv = NA, col = c("grey18", "forestgreen"), cexRow = 0.6,cexCol = 1)
saveWidget(map4, "combined.html")



all<- read.csv("/Users/rckarns8/Qsync/Github/Scripts_Pipelines/87-89_order_profile.csv", header = T, row.names = 1)
map4<-d3heatmap(all, Colv = NA,Rowv = NA, col = c("grey18", "forestgreen"), cexRow = 0.6,cexCol = 1)
saveWidget(map4, "/Users/rckarns8/Qsync/Github/Scripts_Pipelines/combined.html")
