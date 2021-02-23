library(ggtree)
tree <- read.tree(file="/Users/Rachael Karns/Desktop/tree.nwk")
tree2 <- read.tree(file="/Users/Rachael Karns/Downloads/GZXe0N5NdtngoLCUTlZLTw_newick.txt")

plot_tree(tree2)

tr<- tree
d<-as.data.frame(replace_names.tsv)
ggtree(tr) %<+% d + xlim(NA, 5) +
  geom_tiplab(aes(label=paste0(genus), parse=T))

tr2 = rename_taxa(tr, d, label, genus)
write.tree(tr2)
plot_tree(tr2)


t_drop<- c("uncultured", "Calanoida", "Cryomonadida", "Gymnodiniphycidae", "", "Clostridium_sensu_stricto_1", "KD4-96", "Clostridium_sensu_stricto_7", "Clostridium_sensu_stricto_17" )

q5<-drop.tip(tr2, t_drop, trim.internal = TRUE, subtree = FALSE,
             root.edge = 0, rooted = is.rooted(q4), collapse.singles = TRUE,
             interactive = FALSE)	
p<-ggtree(q5)+geom_tiplab()
p1<-ggtree(q5)+geom_text(aes(label=node))
p
p1
ggtree(q5) + geom_point2(aes(subset = sub("/.*", "", label) > 70 | sub(".*/", "", label) > 95))


collapse
library(TreeTools)
test<-CollapseNode(q5,130)
ptest<-ggtree(test)+geom_tiplab()
ptest


d2 = dplyr::mutate(d, newlab = paste(order))
d2
tr3 = rename_taxa(q5, d2, label, newlab)
write.tree(tr3)


