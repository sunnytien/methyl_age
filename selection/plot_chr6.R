source("./selection/gene_network.R")

clusters = get.gene.clusters()

cluster.info = clusters$cluster.info
gene2cluster = clusters$gene2cluster

genes = inner_join(cluster.info, gene2cluster) %>%
  filter(chr=="chr6") %>%
  .$gene

sg = induced.subgraph(g, which(V(g)$name %in% genes))

svg("./figures/chr6_network.svg", 20, 20)
plot(sg,
     vertex.size=1.1,
     vertex.label.cex=0.4,
     vertex.label.color="black",
     vertex.frame.color=NA,
     edge.width=0.5,
     edge.color="black",
     layout=layout.fruchterman.reingold)
dev.off()