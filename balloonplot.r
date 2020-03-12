library(ggpubr)


ggballoonplot(HC_genes, size.range = c(0,4.5), width=500)+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 55),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0))


ggballoonplot(Nitrogen_Genes, size.range = c(0,5), width=500)+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 65),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0))

