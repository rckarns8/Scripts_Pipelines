# Circlepacker package
library(circlepackeR)         


hierarchical_list <- list(name = "World",
                          children = list(
                            list(name = "Halomonas_GOS1",
                                 children = list(
                                   list(name = "Alkane", size = 80),
                                   list(name = "Benzene", size = 64),
                                   list(name = "C1", size = 24),
                                   list(name = "N", size = 52),
                                   list(name = "PAH", size = 78))),
                            list(name = "Alcanivorax_TK-40",
                                 children = list(
                                   list(name = "Alkane", size = 80),
                                   list(name = "Benzene", size = 45),
                                   list(name = "C1", size = 40),
                                   list(name = "N", size = 39),
                                   list(name = "PAH", size = 46))),
                            list(name = "Alcanivorax_TY-5",
                                 children = list(
                                   list(name = "Alkane", size = 100),
                                   list(name = "Benzene", size = 55),
                                   list(name = "C1", size = 36),
                                   list(name = "N", size = 43),
                                   list(name = "PAH", size = 44))),
                            list(name = "Cycloclasticus_TK-8",
                                 children = list(
                                   list(name = "Alkane", size = 28),
                                   list(name = "Benzene", size = 35),
                                   list(name = "C1", size = 80),
                                   list(name = "N", size = 80),
                                   list(name = "PAH", size = 100))),
                            list(name = "Marinobacter_TK-40",
                                 children = list(
                                   list(name = "Alkane", size = 100),
                                   list(name = "Benzene", size = 73),
                                   list(name = "C1", size = 32),
                                   list(name = "N", size = 30),
                                   list(name = "PAH", size = 64))),
                            list(name = "Marinobacter_TT-1",
                                 children = list(
                                   list(name = "Alkane", size = 100),
                                   list(name = "Benzene", size = 91),
                                   list(name = "C1", size = 24),
                                   list(name = "N", size = 22),
                                   list(name = "PAH", size = 74))),
                            list(name = "Pseudoalteromonas_TK-105",
                                 children = list(
                                   list(name = "Alkane", size = 80),
                                   list(name = "Benzene", size = 45),
                                   list(name = "C1", size = 28),
                                   list(name = "N", size = 30),
                                   list(name = "PAH", size = 46))),
                            list(name = "Thalassospira_TK-13(2)",
                                 children = list(
                                   list(name = "Alkane", size = 80),
                                   list(name = "Benzene", size = 64),
                                   list(name = "C1", size = 24),
                                   list(name = "N", size = 57),
                                   list(name = "PAH", size = 64)))
                          ))



circlepackeR(hierarchical_list, color_min = "hsl(202,100%,38%)", 
             color_max = "hsl(128,0%,100%)")








library(ggraph)
library(igraph)
library(tidyverse)
library(viridis)
library(data.tree)

edges <- edges
vertices <- vertices  
mygraph <- graph_from_data_frame(edges, vertices=vertices )


ggraph(mygraph, layout = 'circlepack', weight=vertices$size) + 
  geom_node_circle(aes(fill = as.factor(depth), color = as.factor(depth) )) +
  scale_fill_manual(values=c("0" = "white", "1" = viridis(4)[1], "2" = viridis(4)[2], "3" = viridis(4)[3], "4"=viridis(4)[4])) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black") ) +
  geom_node_text( aes(label=new_label), size=4) +
  theme_void() 




geom.text.size = 4
p<-ggraph(mygraph, layout = 'circlepack', weight=vertices$size) + 
  geom_node_circle(aes(fill = depth)) +
  geom_node_circle() +
  theme_void()+
  geom_node_label(nudge_y = 11.75, aes(label=vertices$new_label, size=18), size=geom.text.size)+
  scale_fill_distiller(palette = "GnBu")
p



