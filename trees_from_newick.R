#Load Libraries

library(ggtree)
library(TreeTools)
library(treeio)
library(phytools)
library(ggplot2)



#Load in newick tree
tree <- read.tree(file="/Users/Rachael Karns/Desktop/BOEM/Sed4/sed4/tree.nwk")

#Load in a file containing the following headers: 
#label  domain  phylum class order family genus species
#It is easy to make this from the qiime output in excel- separate by ';' and remove *__
#Load/import this file in
d<-as.data.frame(replace_names.tsv)

#Plot the tree as it is

ggtree(tree) %<+% d + xlim(NA, 5) +
  geom_tiplab(aes(label=paste0(family), parse=T))


d2 = dplyr::mutate(d, newlab = paste(family,genus, sep='_'))
d2

d2 = dplyr::mutate(d, newlab = paste(genus))

tr3 = rename_taxa(tree, d2, label, newlab)
write.tree(tr3)





## here are our family names

tips<-tr3$tip.label
family<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
## here are our families
family

ii<-sapply(family,function(x,y) grep(x,y)[1],y=tips)
treetest<-drop.tip(tr3,setdiff(tr3$tip.label,tips[ii]))



t_drop<- c("uncultured", "Bacillus","Calanoida", "Cryomonadida", "Gymnodiniphycidae", "", "Clostridium_sensu_stricto_1", "KD4-96", "Clostridium_sensu_stricto_7", "Clostridium_sensu_stricto_17" )

tree_prune<-drop.tip(treetest, t_drop, trim.internal = TRUE, subtree = FALSE,
                     root.edge = 0, rooted = is.rooted(treetest), collapse.singles = TRUE,
                     interactive = FALSE)	
fig.width=7
ggtree(tree_prune, layout="circular", open.angle=10)+
  geom_tiplab(size=2.5,hjust = -.5, align = FALSE, linetype = "dotted", linesize = .5)+

  
#Save the figure as a pdf
ggsave("/Users/Rachael Karns/Desktop/test.pdf", width = 50, height = 50, units = "cm", limitsize = FALSE)

  



