#Heatmap generator
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: December 3, 2019 by Rachael Storo

#Purpose: Create a script for visualizing relative abundance tables as a heat map


library(d3heatmap)
library(htmlwidgets)

all<- read.csv("/work/sbjlab/rck/Sed_Assemblies/combined_87-89-91_Metaphlan_profile.csv", header = T, row.names = 1)
map4<-d3heatmap(all, Colv = NA,Rowv = NA, col = c("grey18", "forestgreen"), cexRow = 0.6,cexCol = 1, xaxis_height = 100, yaxis_width = 275 )
saveWidget(map4, "combined.html")