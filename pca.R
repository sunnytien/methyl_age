library("dplyr")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "test4")

sites = beta %>%
  group_by(gsm.id) %>%
  summarize(n=n()) 

beta.subsample = beta %>%
  semi_join(sites)

