library("dplyr")
library("tidyr")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

sites = beta %>%
  select(Probe) %>%
  collect() %>%
  sample_n(1000)

x = beta %>%
  semi_join(sies, copy=T) %>%
  collect %>%
  as.matrix

p = pca(x, scale="none", center=T)

save(p, file="./data/pca.Rdata")
