## script to count the number of probes present in BMIQ files

library("dplyr")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "test4")

probes.present = beta %>%
  group_by(gsm.id) %>%
  summarize(n.probes=n()) %>%
  collect

save(probes.present, file="probes.present.Rdata")