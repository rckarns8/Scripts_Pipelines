#Make a heatmap from MEBS data 

#Import data
library(tidyverse)
library("RColorBrewer")

samp.with.rownames <- data.frame(Oil_Fluff_FDR0.001[,-1], row.names=Oil_Fluff_FDR0.001[,1])
samp.with.rownames
as.matrix(samp.with.rownames)

par(oma=c(3,3,3,3)) # all sides have 3 lines of space
par(mar=c(5,4,4,2) + 0.1)
heatmap.2(as.matrix(samp.with.rownames[,-1]), lhei=c(1,4.5),lwid=c(1,3),density.info="none",dendrogram='none',Colv=FALSE, Rowv= FALSE,key=TRUE, 
            symkey=FALSE,tracecol=NA,scale="column",keysize=1, na.color = "white",col=brewer.pal(n=11,name="PRGn"))
