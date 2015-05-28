library("tidyr")
library("dplyr")
library("data.table")


# this script processes the results from the lmers 
# and combined everything into a single data.frame
# for easy analysis 

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

# ## process anova results
# 
# anova.files = list.files("./data/anovas", full.names=T)
# 
# anova.results = anova.files %>%
#   lapply(process.anova) 
# 
# ncols = sapply(anova.results, ncol)
# 
# anova.result = anova.results[ncols == 8] %>%
#   rbindlist
# 

## process lmer results

summary.files = list.files("./data/models", full.names=T)

summary.results = summary.files %>%
  lapply(process.summary) 

ncols = sapply(summary.results, ncol) 
summary.result = summary.results[ncols==8] %>% # if ncol != 8, something went wrong with model fitting
  rbindlist

# effects of aging
age = summary.result %>%
  filter(variable=="age") %>%
  arrange(Estimate) %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr"))

# effect of AFR vs EUR
afr = summary.result %>%
  filter(variable=="predicted.ancestryEUR") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  arrange(Estimate)

# difference between solid and liquid tissues
tissue = summary.result %>%
  filter(variable=="tissue_state1") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  arrange(Estimate)

# interactions...
age.pop = summary.result %>%
  filter(variable=="age:predicted.ancestryEUR") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  arrange(Estimate)

tissue.pop = summary.result %>%
  filter(variable=="tissue_state1:predicted.ancestryEUR") %>%
  mutate(q.value=p.adjust(`Pr(>|t|)`, method="fdr")) %>%
  arrange(Estimate)

save(age, file="./data/age.Rdata")
save(tissue, file="./data/tissue.Rdata")
save(afr, file="./data/afr.Rdata")