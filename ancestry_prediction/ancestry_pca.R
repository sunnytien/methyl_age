library("dplyr")
library("randomForest")
library("tidyr")
library("pcaMethods")

load("./data/sample.info.Rdata")
load("./data/rf.Rdata")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

probes = rownames(importance(rf))

b = beta %>%
  filter(Probe %in% probes) %>%
  collect

b.matrix = b %>%
  select(-Probe) %>%
  as.matrix %>%
  t 

good.samples = sample.info %>% 
  filter(!(tissue %in% c("Lymphocytes", "Lymphoblasts"))) %>%
  filter(gsm.ind %in% rownames(b.matrix))

p = pca(b.matrix[good.samples$gsm.id,1:100], center=T, scale="none")
