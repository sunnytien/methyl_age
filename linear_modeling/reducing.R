library("data.table")

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

save(anova.results, file="./Rmd/gene_selection/anova.results.Rdata")

## process lmer results

summary.files = list.files("./data/models", full.names=T)

summary.results = files %>%
  lapply(process.summary)

ncols = sapply(summary.results, ncol)

summary.result = summary.results[ncols==7] %>%
  rbindlist

save(lmer.results, file="./Rmd/gene_selection/lmer.results.Rdata")

age = summary.result %>%
  filter(variable=="age.normed")

age.est = age$Estimate
names(age.est) = age$gene

pop.est = pop$Estimate
names(pop.est) = pop$gene
