#Heatmap generator
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: December 3, 2019 by Rachael Storo

#Purpose: Create a script for visualizing relative abundance tables as a heat map


library(d3heatmap)
hmm3<- read.csv("/work/sbjlab/rck/Sed_Assemblies/AT26-13-87/87_order_profile.csv", header = T, row.names = 1)
d3heatmap(hmm3, Colv = NA,Rowv = NA, col = c("grey28", "darkseagreen"), cexRow = 0.6,cexCol = 1)

hmm3<- read.csv("/work/sbjlab/rck/Sed_Assemblies/AT26-13-89/89_order_profile.csv", header = T, row.names = 1)
d3heatmap(hmm3, Colv = NA,Rowv = NA, col = c("grey28", "darkseagreen"), cexRow = 0.6,cexCol = 1)

hmm3<- read.csv("/work/sbjlab/rck/Sed_Assemblies/AT26-13-91/91_order_profile.csv", header = T, row.names = 1)
d3heatmap(hmm3, Colv = NA,Rowv = NA, col = c("grey28", "darkseagreen"), cexRow = 0.6,cexCol = 1)


