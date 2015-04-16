library("GenomicRanges")
library("FDb.InfiniumMethylation.hg19")

source("./selection/gene_network.R")

get.d.statistic = function(scores){
  ks = ks.test(scores, sample(null, 1e4))
  ks$statistics
}

load("./data/afr.Rdata")
load("./data/xpehh.Rdata")

gene.granges = get.gene.granges(afr$gene)

hits = findOverlaps(xpehh, gene.granges)

pairs = data.frame(snpID=xpehh$mcols.snpID[queryHits(hits)],
                   scores=xpehh$mcols.score[queryHits(hits)],
                   gene=gene.granges$mcols.gene[subjectHits(hits)],
                   stringsAsFactors=F)

null = xpehh$mcols.snpID[queryHits(hits) %>% unique]

d.statistics = pairs %>%
  group_by(gene) %>%
  do(d=get.d.statistic(.$scores))