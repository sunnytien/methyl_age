library("dplyr")
library("tidyr")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "test4")

sites = beta %>%
  group_by(probe.id) %>%
  summarize(n=n()) %>%
  collect

sites.sample = sites %>%
  filter(n == max(n)) %>%
  sample_n(10) %>%
  select(probe.id)

beta.subsample = beta %>%
  semi_join(sites.sample, copy=T) %>%
  collect %>%
  spread(probe.id, beta) 

x = beta.subsample %>%
  select(starts_with("cg")) %>%
  as.matrix
rownames(x) = beta.subsample$gsm.id

p = pca(x, scale="none", center=T)

save(p, file="./data/pca.Rdata")
