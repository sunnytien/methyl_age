source("./selection/gene_network.R")

load("./data/afr.Rdata")
load("./data/xpehh.Rdata")

gene.granges = get.gene.granges(afr$gene)

hits = findOverlaps(xpehh, gene.granges)

pairs = data.frame(snpID=xpehh$mcols.snpID[queryHits(hits)],
                   scores=xpehh$mcols.score[queryHits(hits)],
                   gene=gene.granges$mcols.gene[subjectHits(hits)],
                   stringsAsFactors=F)

probs = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999)

average.xpehh = pairs %>%
  mutate(scores=abs(scores)) %>%
  group_by(gene) %>%
  summarize(mean=mean(scores),
            min=min(scores),
            max=max(scores),
            sd=sd(scores),
            q=list(quantile(scores, probs=probs)))


tmp = inner_join(average.xpehh,
                 afr)

z2 = sapply(tmp$q, function(x) x[2])
z8 = sapply(tmp$q, function(x) x[8])

z = ifelse(z2 > z8, z2, z8)


pos = tmp %>% 
  filter(mean > 0)

neg = tmp %>%
  filter(mean < 0)

### fst

hits2 = findOverlaps(fst, gene.granges)

pairs = data.frame(snpID=fst$mcols.snpID[queryHits(hits2)],
                   scores=fst$mcols.score[queryHits(hits2)],
                   gene=gene.granges$mcols.gene[subjectHits(hits2)],
                   stringsAsFactors=F)

probs = c(0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999)

average.fst = pairs %>%
  group_by(gene) %>%
  summarize(mean=mean(scores),
            min=min(scores),
            max=max(scores),
            sd=sd(scores),
            q=list(quantile(scores, probs=probs)))

tmp = inner_join(average.fst,
                 afr)

z1 = sapply(tmp$q, function(x) x[1])
z2 = sapply(tmp$q, function(x) x[2])
z3 = sapply(tmp$q, function(x) x[3])

z8 = sapply(tmp$q, function(x) x[8])
