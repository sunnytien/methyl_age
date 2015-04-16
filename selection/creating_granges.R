library("data.table")
library("dplyr")
library("GenomicRanges")

fst.dt = fread("./data/FstGLOB_CEU_u_YRI_u_CHB.whole_genome.pvalues")

fst = GRanges(seqnames=paste("chr", fst.dt$chromosome, sep=""),
              ranges=IRanges(fst.dt$position,
                             fst.dt$position),
              strand="*",
              mcols=DataFrame(snpID=fst.dt$snpID,
                              score=fst.dt$score,
                              pvalue=fst.dt$pvalue))

save(fst, file="./data/fst.Rdata")

xpehh.dt = fread("./data/XPEHH_CEU_vs_YRI.whole_genome.pvalues")

xpehh = GRanges(seqnames=paste("chr", xpehh.dt$chromosome, sep=""),
                ranges=IRanges(xpehh.dt$position,
                               xpehh.dt$position),
                strand="*",
                mcols=DataFrame(snpID=xpehh.dt$snpID,
                                score=xpehh.dt$score,
                                pvalue=xpehh.dt$pvalue))

save(xpehh, file="./data/xpehh.Rdata")