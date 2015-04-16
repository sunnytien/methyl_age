library("data.table")
library('dplyr')

source("./selection/gene_network.R")

load("./data/afr.Rdata")

fst = fread("./data/FstGLOB_CEU_u_YRI_u_CHB.whole_genome.pvalues")

fst = GRanges(seqnames=paste("chr", fst$chromosome),
              ranges=IRanges(fst$position,
                              fst$position),
              strand=*,
              mcols=DataFrame(score=fst$score,
                              pvalue=fst$pvalue))


## xpehh = fread("./data")



## finding snps that intersect with all promoters
all = get.gene.clusters()



regions = GRanges()

## intersect SNPs with top
top = get.gene.clusters()
