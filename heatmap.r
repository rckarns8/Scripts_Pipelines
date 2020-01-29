#Heatmap generator
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: December 3, 2019 by Rachael Storo

#Purpose: Create a script for visualizing relative abundance tables as a heat map


library(d3heatmap)
library(htmlwidgets)

all<- read.csv("~/Desktop/Analysis/Projects/Comparative_Genomics_Papers/Gulf/Other_Genes.csv", header = T, row.names = 1)
map4<-d3heatmap(all, Colv = NA,Rowv = NA, col = c("grey34", "darkseagreen"), xaxis_height = 200, yaxis_width = 190, xaxis_font_size= 16, yaxis_font_size= 12 )
map4
saveWidget(map4, "combined.html")
