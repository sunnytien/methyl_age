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

process.lmer = function(x){
  
  gene = gsub("./data/models/", "", x) %>%
    gsub(".Rdata", "", .)
  
  m = get(load(x))
  
  df = data.frame(gene=gene,
                  variable=names(fixef(m)),
                  coef=fixef(m))
}

library("lme4")
library("lmerTest")
library("dplyr")

## process anova results

files = list.files("./data/anovas", full.names=T)

anova.results = files %>%
  lapply(process.anova) 

ncols = sapply(anova.results, ncol)

anova.results = anova.results[ncols==8] %>%
  do.call(rbind, .)
  
save(anova.results, file="./data/anova.results.Rdata")

## process lmer results

files = list.files("./data/models", full.names=T) 

lmer.results = files %>%
  lapply(process.lmer) %>%
  do.call(rbind, .)