library("GenomicRanges")
library("FDb.InfiniumMethylation.hg19")
library("ggplot2")

source("./selection/gene_network.R")

get.d.statistic = function(scores){
  ks = ks.test(scores, sample(null, 1e4))
  ks$statistic
}

load("./data/afr.Rdata")
load("./data/xpehh.Rdata")
load("./data/fst.Rdata")

gene.granges = get.gene.granges(afr$gene)

hits = findOverlaps(xpehh, gene.granges)

pairs = data.frame(snpID=xpehh$mcols.snpID[queryHits(hits)],
                   scores=xpehh$mcols.score[queryHits(hits)],
                   gene=gene.granges$mcols.gene[subjectHits(hits)],
                   stringsAsFactors=F)

null = xpehh$mcols.score[hits %>% queryHits %>% unique]

summary = pairs %>%
  group_by(gene) %>%
  do(d=get.d.statistic(.$scores))

d = unlist(summary$d)

summary$D = d

tmp = inner_join(summary, afr)

ggplot(tmp, aes(D, Estimate)) +
  geom_point() +
  geom_smooth(method="lm")

plot(tmp$D, abs(tmp$Estimate),
    pch=20,
    bty="n",
    xlab="Difference from Expected XP-EHH Distribution",
    ylab="Methyatlion difference between populations")

m = lm(abs(Estimate) ~ D, data=tmp)
abline(m, col="red")

## outlier proportion
cutoff = quantile(null, probs=c(0.025, 0.975))

outliers = pairs %>%
  group_by(gene) %>%
  summarize(low=sum(scores < cutoff[1])/n(),
            high=sum(scores > cutoff[2])/n()) %>%
  mutate(total=low+high) %>%
  inner_join(afr)

plot(outliers$total, abs(outliers$Estimate))
cor(outliers$low, outliers$Estimate)

## fst proportion

hits = findOverlaps(fst, gene.granges)

pairs = data.frame(snpID=fst$mcols.snpID[queryHits(hits)],
                   scores=fst$mcols.score[queryHits(hits)],
                   gene=gene.granges$mcols.gene[subjectHits(hits)],
                   stringsAsFactors=F)

null = fst$mcols.score[hits %>% queryHits %>% unique]
cutoff = quantile(null, prob=0.99)

outliers = pairs %>%
  group_by(gene) %>%
  summarize(p=mean(scores)) %>%
  inner_join(afr)