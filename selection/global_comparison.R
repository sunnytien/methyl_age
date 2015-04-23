library("FDb.InfiniumMethylation.hg19")
library("data.table")
library('dplyr')

source("./selection/gene_network.R")

select = dplyr::select

load("./data/afr.Rdata")
load("./data/fst.Rdata")
load("./data/xpehh.Rdata")

## finding snps that intersect with all promoters
all = get.gene.clusters()

all.xpehh.hits = findOverlaps(xpehh, all$grange)
all.xpehh = xpehh[queryHits(all.xpehh.hits)] %>%
  as.data.frame

all.fst.hits = findOverlaps(fst, all$grange)
all.fst = fst[queryHits(all.fst.hits)] %>%
  as.data.frame

## intersect SNPs with top
top.genes = afr %>%
  arrange(-abs(`t value`)) %>%
  .$gene %>%
  head(50)

top = get.gene.clusters(top.genes)

top.xpehh.hits = findOverlaps(xpehh, top$grange)
top.xpehh = xpehh[queryHits(top.xpehh.hits)] %>%
  as.data.frame

top.fst.hits = findOverlaps(fst, top$grange)
top.fst = fst[queryHits(top.fst.hits)] %>%
  as.data.frame

hist(all.xpehh$mcols.score)
hist(top.xpehh$mcols.score)

pdf("./figures/selection.pdf")

boxplot(top.xpehh$mcols.score, all.xpehh$mcols.score, outline=F,
        names=c("Top DM Promoters", "All Promoters"),
        frame=F,
        ylab="XP-EHH between CEU and YRI")



boxplot(top.fst$mcols.score, all.fst$mcols.score, outline=F,
        names=c("Top DM Promoters", "All Promoters"),
        frame=F,
        ylab="Fst between CEU and YRI")

dev.off()

var.test(top.xpehh$mcols.score, all.xpehh$mcols.score)
wilcox.test(top.fst$mcols.score, all.fst$mcols.score)

