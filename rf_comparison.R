library("randomForest")
library("dplyr")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

mm = read.table("./data/probes_MM.txt")

old.rf = get(load("./data/rf.Rdata"))
new.rf = get(load("./data/ancestry_rf.Rdata"))

old.imp = importance(old.rf) %>%
  as.data.frame %>%
  mutate(probe.id=rownames(.)) %>%
  select(probe.id, old.gini=MeanDecreaseGini)

new.imp = importance(new.rf) %>%
  as.data.frame %>%
  mutate(probe.id=rownames(.)) %>%
  select(probe.id, new.gini=MeanDecreaseGini)

imp = full_join(old.imp, new.imp) %>%
  mutate(diff = old.gini - new.gini) %>%
  arrange(diff)

imp2 = imp %>%
  filter(is.na(diff)) %>%
  arrange(-new.gini)


middle.rf = randomForest(x.training[,rownames(importance(old.rf))], y.training)
