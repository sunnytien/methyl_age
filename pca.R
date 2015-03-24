library("dplyr")
library("tidyr")
library("pcaMethods")

load("./data/sample.info.Rdata")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

sites = beta %>%
  select(Probe) %>%
  collect() %>%
  sample_n(1000) %>%
  .$Probe

x = beta %>%
  filter(Probe %in% sites) %>%
  select(starts_with("GSM")) %>%
  collect %>%
  as.matrix %>%
  t
  

no_outliers = sample.info %>%
  filter(tissue != "Lymphocytes")

blood = sample.info %>%
  filter(tissue == "Whole Blood")

brain = sample.info %>%
  filter(tissue == "Brain")

p.all = pca(x, scale="none", center=T)
p.nooutliers = pca(x[rownames(x) %in% no_outliers$gsm.id, ], scale="none", center=T)
p.blood = pca(x[rownames(x) %in% blood$gsm.id, ], scale="none", center=T)
p.brain = pca(x[rownames(x) %in% brain$gsm.id, ], scale="none", center=T)

save(p.all,
     p.nooutliers,
     p.blood,
     p.brain, file="./data/pca.Rdata")
