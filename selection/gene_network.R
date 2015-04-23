library("FDb.InfiniumMethylation.hg19")
library("dplyr")
library("igraph")
library("RColorBrewer")
library("GenomicRanges")

get.gene.clusters = function(genes=NULL, d=40e3){
  
  load("./data/transcripts.Rdata")
  
  grl <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
  gr <- unique(resize(unlist(grl), 1, fix = "start"))
  
  names(gr) = NULL
  
  tss = as.data.frame(gr) %>%
    select(chr=seqnames, tss=start, tx=tx_name) %>%
    inner_join(transcripts)
  
  tss = if(!is.null(genes)) tss %>% filter(gene %in% genes) else
    tss
  
  pairs = inner_join(tss %>% select(chr=chr,
                                    tss1=tss,
                                    tx1=tx,
                                    gene1=gene),
                     tss %>% select(chr=chr,
                                    tss2=tss,
                                    tx2=tx,
                                    gene2=gene)) %>%
    filter(tss1 > tss2) %>%
    filter(abs(tss1 - tss2) < 2*d)
  
  g = graph.empty(directed=F) + 
    vertices(tss$gene) +
    edges(cbind(pairs$gene1, pairs$gene2) %>% t %>% as.vector)
  
  cl = clusters(g)
  
  ## cluster membership
  gene2cluster = data.frame(gene=V(g)$name,
                            cluster=cl$membership,
                            stringsAsFactors=F)
  
  cluster.info = gene2cluster %>%
    inner_join(tss) %>%
    group_by(cluster) %>%
    summarize(min=min(tss), max=max(tss), chr=chr[1], n=n()) %>%
    mutate(start=min-d, end=max+d)
  
  V(g)$color = brewer.pal(12, "Set3")[cl$membership %% 12 + 1]
  
  granges = GRanges(seqnames=cluster.info$chr,
                   ranges=IRanges(cluster.info$start,
                                  cluster.info$end),
                   strand="*",
                   mcols=DataFrame(cluster=cluster.info$cluster))
  
  list(g=g,
       cluster.info=cluster.info,
       gene2cluster=gene2cluster,
       grange=granges)
}

get.gene.granges = function(genes=NULL, d=40e3){
  
  load("./data/transcripts.Rdata")
  
  grl <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
  gr <- unique(resize(unlist(grl), 1, fix = "start"))
  
  names(gr) = NULL
  
  tss = as.data.frame(gr) %>%
    select(chr=seqnames, tss=start, tx=tx_name) %>%
    inner_join(transcripts)
  
  tss = if(!is.null(genes)) tss %>% filter(gene %in% genes) else
    tss
  
  GRanges(seqnames=tss$chr,
          ranges=IRanges(tss$tss-d,
                         tss$tss+d),
          strand="*",
          mcols=DataFrame(gene=tss$gene))
  
}
