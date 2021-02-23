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




q<-drop.tip(tr2, "", trim.internal = TRUE, subtree = FALSE,
            root.edge = 0, rooted = is.rooted(tr2), collapse.singles = TRUE,
            interactive = FALSE)
q
q2<-drop.tip(q, "uncultured", trim.internal = TRUE, subtree = FALSE,
            root.edge = 0, rooted = is.rooted(q), collapse.singles = TRUE,
            interactive = FALSE)
q3<-drop.tip(q2, "Calanoida", trim.internal = TRUE, subtree = FALSE,
             root.edge = 0, rooted = is.rooted(q2), collapse.singles = TRUE,
             interactive = FALSE)
q4<-drop.tip(q3, "Cryomonadida", trim.internal = TRUE, subtree = FALSE,
             root.edge = 0, rooted = is.rooted(q3), collapse.singles = TRUE,
             interactive = FALSE)
q5<-drop.tip(q4, "Gymnodiniphycidae", trim.internal = TRUE, subtree = FALSE,
             root.edge = 0, rooted = is.rooted(q4), collapse.singles = TRUE,
             interactive = FALSE)	


	
d2 = dplyr::mutate(d, newlab = paste(order))
d2
tr3 = rename_taxa(q5, d2, label, newlab)
write.tree(tr3)
p<-ggtree(q5)+geom_tiplab()
p
ggtree(tr3) + geom_point2(aes(subset = sub("/.*", "", label) > 70 | sub(".*/", "", label) > 95))
