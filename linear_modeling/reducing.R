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
    mutate(variable=rownames(.)) %>%
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

summary.result = summary.results[ncols==7] %>%
  rbindlist

age = summary.result %>%
  filter(variable=="age.normed") %>%
  arrange(Estimate)
save(age, file="./data/age.Rdata")

age.est = age$Estimate
names(age.est) = age$gene

age.gsea = gsea(age.est, gmt="ReactomePathways.gmt")
save(age.gsea, file="./data/age.gsea.Rdata")


tissue = summary.result %>%
  filter(variable=="tissue_state1")

tissue.est = tissue$Estimate
names(tissue.est) = tissue$gene

tissue.gsea = gsea(tissue.est, gmt="ReactomePathways.gmt")
save(tissue.gsea, file="./data/tissue.gsea.Rdata")


pop = anova.result %>%
  filter(variable=="predicted.ancestry")

pop.est = pop$`F.value`
names(pop.est) = pop$gene
pop.gsea = gsea(pop.est, gmt="ReactomePathways.gmt")

asn = summary.result %>%
  filter(variable=="predicted.ancestryASN")
asn.est = asn$Estimate
names(asn.est) = asn$gene
asn.gsea = gsea(asn.est, gmt="ReactomePathways.gmt")

afr = summary.result %>%
  filter(variable=="predicted.ancestryEUR")
afr.est = afr$Estimate
names(afr.est) = afr$gene
afr.gsea = gsea(afr.est, gmt="ReactomePathways.gmt", output.dir=".")
