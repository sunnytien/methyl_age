library("dplyr")
library("tidyr")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "test4")

sites = beta %>%
  select(probe.id) %>%
  distinct %>%
  head(1000)

sites.sample = sites %>%
  sample_n(1000) 

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
