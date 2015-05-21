library("tidyr")
library("dplyr")
library("data.table")

source("~/Projects/GSEA/R/gsea.R")

### this script loads the results of the anova 
### and turns it into a single data.frame for
### graphical display by Rmd files

process.anova = function(x){
  
  gene = gsub("./data/anovas/", "", x) %>%
    gsub(".Rdata", "", .)
  
  df = get(load(x))
  
  df %<>% 
    mutate(variable=rownames(.)) %>%
    mutate(gene=gene)

  return(df)
}

process.summary = function(x){
  
  gene = gsub("./data/models/", "", x) %>%
    gsub(".Rdata", "", .)
  
  m = get(load(x))
  
  df = m %>% 
    as.data.frame %>%
    mutate(gene=gene)
  return(df)
}

library("lme4")
library("lmerTest")
library("dplyr")

## process anova results

anova.files = list.files("./data/anovas", full.names=T)

anova.results = anova.files %>%
  lapply(process.anova) 

ncols = sapply(anova.results, ncol)

anova.result = anova.results[ncols == 8] %>%
  rbindlist


## process lmer results

summary.files = list.files("./data/models", full.names=T)

summary.results = summary.files %>%
  lapply(process.summary) 

ncols = sapply(summary.results, ncol) 
summary.result = summary.results[ncols==8] %>%
  rbindlist

age = summary.result %>%
  filter(variable=="age") %>%
  arrange(Estimate) %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr"))

age %>%
  filter(q.value < 0.01) %>%
  select(gene) %>%
  write.table(quote=F, row.names=F, col.names=F, sep=",",
              file="~/Desktop/age_genes.csv")


save(age, file="./data/age.Rdata")

age.est = age$Estimate
names(age.est) = age$gene

age.gsea = gsea(age.est, gmt="ReactomePathways.gmt")
save(age.gsea, file="./data/age.gsea.Rdata")


age.genes = age %>%
  filter(q.value < 0.01) 

write.table(age.genes %>% select(gene, dBeta) %>% filter(dBeta > 0),
            sep="\t",
            row.names=F,
            col.names=F,
            file="~/Desktop/age_up.txt")

write.table(age.genes %>% select(gene, dBeta) %>% filter(dBeta < 0),
            sep="\t",
            row.names=F,
            col.names=F,
            file="~/Desktop/age_down.txt")



afr = summary.result %>%
  filter(variable=="predicted.ancestryEUR") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr"))

afr.est = afr$Estimate
names(afr.est) = afr$gene
afr.gsea = gsea(afr.est, gmt="ReactomePathways.gmt", output.dir=".")


age.pop = summary.result %>%
  filter(variable=="age:predicted.ancestryEUR") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  arrange(Estimate)

tissue.pop = summary.result %>%
  filter(variable=="tissue_state1:predicted.ancestryEUR") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  arrange(Estimate)

tissue.pop.est = tissue.pop$Estimate
names(tissue.pop.est) = tissue.pop$gene

tissue.pop.gsea = gsea(tissue.pop.est, gmt="ReactomePathways.gmt")

tissue = summary.result %>%
  filter(variable=="tissue_state1") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  arrange(Estimate)

get.sig = function(x) x %>%
  filter(q.value < 0.05) %>%
  filter(abs(dBeta) > 0.05) %>%
  .$gene %>%
  unique

tissue.genes = get.sig(tissue)
afr.genes = get.sig(afr)
age.genes = get.sig(age)

pdf("~/Desktop/venn.pdf")
draw.triple.venn(length(tissue.genes),
                 length(age.genes),
                 length(afr.genes),
                 length(intersect(tissue.genes, age.genes)),
                 length(intersect(age.genes, afr.genes)),
                 length(intersect(tissue.genes, afr.genes)),
                 length(intersect(intersect(tissue.genes, age.genes), afr.genes)),
                 col=brewer.pal(3, "Set3"),
                 fill=brewer.pal(3, "Set3"),
                 euler.d=T,
                 scaled=T,
                 category=c("Tissue DMPs", "Age DMPs", "Ancestry DMPs"),
                 cat.cex=1.3)
dev.off()

load("./data/model_data/ACKR1.Rdata")

library("ggplot2")
cg = data %>% filter(Probe=="cg04922029")

ggplot(cg, aes(predicted.ancestry, M, fill=predicted.ancestry)) + 
  geom_boxplot(outlier.colour=NA) + 
  facet_grid(Probe ~ tissue_state) +
  theme_bw() +
  ylim(-2.5, 3) +
  xlab("") +
  ylab("M value") +
  ggtitle("DARC CpG Site") +
  guides(color=F) +
  theme(plot.title=element_text(size=rel(2))) + 
  theme(axis.title.y=element_text(size=rel(1.3)))

