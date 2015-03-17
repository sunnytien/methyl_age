library("pcaMethods")
library("dplyr")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "test4")

sites = beta %>%
  group_by(probe.id) %>%
  summarize(n=n()) %>%
  ungroup %>%
  filter(n==max(n)) %>%
  sample_n(1000)

beta.subsample = beta %>%
  semi_join(sites)

