library("caret")
library("randomForest")
library("glmnet")
library("gbm")
library("kernlab")
library("dplyr")
library("doParallel")

age.trans = function(x) sapply(x, function(y) if(y < 20) log((y + 1)/21) + 1 else (y+1)/21)

load("./data/sample.info.Rdata")
load("./data/probe.info.Rdata")

### DATA PREPARATION ###

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

grm2 = probe.info %>%
  filter(nearestGeneSymbol=="GRM2")

b.tmp = beta %>%
  filter(Probe %in% grm2$Probe) %>%
  collect %>%
  gather(gsm.id, beta, starts_with("GSM")) %>% # transposing
  spread(Probe, beta) # this may be more efficient with matrix stuff

data = sample.info %>%
  select(gsm.id, age) %>%
  filter(!is.na(age)) %>%
  mutate(age.normed=age.trans(age)) %>%
  select(-age) %>%
  inner_join(b.tmp) %>%
  sample_n(500)

x = data %>%
  select(starts_with("cg")) %>%
  as.matrix

y = data$age.normed

### MODEL RUNNING ###

registerDoParallel(5)

ctrl = trainControl(method="repeatedcv",
                    number=3,
                    repeats=1,
                    savePredictions=T,
                    allowParallel=T)

m = train(x, y,
          trControl=ctrl)
