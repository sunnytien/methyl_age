library("dplyr")
library("randomForest")
library("tidyr")
library("pcaMethods")
library("ggvis")

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
colnames(b.matrix) = b$Probe

m.matrix = log(b.matrix / (1 - b.matrix))

good.samples = sample.info %>% 
  filter(!(tissue %in% c("Lymphocytes", "Lymphoblasts"))) %>%
  filter(gsm.id %in% rownames(b.matrix))

m.p = pca(m.matrix[good.samples$gsm.id, sample(1:ncol(b.matrix), 1000)], 
          center=T, 
          scale="none",
          nPcs=4)

b.p = pca(b.matrix[good.samples$gsm.id, sample(1:ncol(b.matrix), 1000)], 
          center=T, 
          scale="none",
          nPcs=4)


m.p.scores = m.p %>%
  scores %>%
  as.data.frame %>%
  mutate(gsm.id=rownames(.)) %>%
  inner_join(sample.info)

b.p.scores = b.p %>%
  scores %>%
  as.data.frame %>%
  mutate(gsm.id=rownames(.)) %>%
  inner_join(sample.info)

imp = rf %>%
  importance %>%
  as.data.frame %>%
  mutate(gsm.id=rownames(.)) %>%
  arrange(-MeanDecreaseGini) %>%
  head(100) %>%
  .$gsm.id

b.p.imp = pca(b.matrix[good.samples$gsm.id, imp], 
          center=T, 
          scale="none",
          nPcs=4)

b.p.imp.scores = b.p.imp %>%
  scores %>%
  as.data.frame %>%
  mutate(gsm.id=rownames(.)) %>%
  inner_join(sample.info) %>%
  inner_join(predicted.ancestry %>%
               select(gsm.id, predicted.ancestry))

plot(b.p.imp.scores$PC1,
     col=factor(b.p.imp.scores$ancestry))

pdf("./figures/ancestry.pdf")

b.p.scores %>%
  ggvis(~PC1, ~PC2, fill:="grey",opacity:=0.2) %>%
  layer_points() %>%
  add_data(b.p.scores %>%
             filter(!is.na(ancestry))) %>%
  layer_points(fill=~ancestry, opacity:=1)

b.p.imp.scores %>%
  ggvis(~PC1, ~PC2, fill:="grey",opacity:=0.2) %>%
  layer_points() %>%
  add_data(b.p.imp.scores %>%
             filter(!is.na(ancestry))) %>%
  layer_points(fill=~ancestry, opacity:=1)

b.p.imp.scores %>%
  ggvis(~PC1, ~PC2, fill:="grey", opacity:=0.2) %>%
  layer_points() %>%
  add_data(b.p.imp.scores %>%
             filter(!is.na(ancestry)) %>%
             filter(!is.na(age)) %>%
             filter(age < 18)) %>%
  layer_points(fill=~ancestry, opacity:=1)

b.p.imp.scores %>%
  ggvis(~PC1, ~PC2, fill:="grey", opacity:=0.2) %>%
  layer_points() %>%
  add_data(b.p.imp.scores %>%
             filter(!is.na(ancestry)) %>%
             filter(!is.na(age)) %>%
             filter(age < 18)) %>%
  layer_points(fill=~predicted.ancestry, opacity:=1)

dev.off()

