library("dplyr")
library("tidyr")
library("ggvis")
library("randomForest")
library("magrittr")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450k.db")

load("./data/sample.info.Rdata")

## preparing labels for training set

get.ancestry = function(x){
  if(grepl("African", x)) return("AFR")
  if(grepl("Caucasian", x)) return("EUR")
  if(grepl("Han|Japanese",x )) return("ASN")
  return(NA)
}

training.ancestry = read.delim("./data/training_ancestry.txt",
                               header=T,
                               stringsAsFactors=F) %>%
  mutate(ancestry=sapply(ancestry, get.ancestry)) %>%
  mutate(gsm.id=gsub(" ", "", gsm.id))

load("./data/aimms.Rdata")
names(aimms) = NULL

## preparing beta-values

chrXY = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations %>%
  as.data.frame %>%
  mutate(Probe=rownames(.)) %>%
  filter(chr %in% c("chrX", "chrY"))

multiple_match = read.table("./data/probes_MM.txt") %>%
  select(Probe=V1)

IlluminaHumanMethylation450kDMR %<>% as.data.frame

aimms = aimms[!(aimms %in% chrXY$Probe)]
aimms = aimms[!(aimms %in% multiple_match$Probe)]
aimms = aimms[!(aimms %in% IlluminaHumanMethylation450kDMR$Probe_ID)]

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

type = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest %>%
  as.data.frame %>%
  mutate(Probe=rownames(.)) %>%
  mutate(type_rg=paste(Type, Color, sep="_")) %>%
  select(Probe, type_rg) %>%
  filter(Probe %in% aimms)

b.wide = beta %>%
  filter(Probe %in% aimms) %>%
  inner_join(type, copy=T) %>%
  collect 

b.thin = b.wide %>%
  gather(gsm.id, beta, -Probe, -type_rg)

## quantile normalizing beta-values

b.thin.normed = b.thin %>%
  group_by(gsm.id) %>%
  mutate(q=rank(beta)) %>%
  ungroup

b.wide.normed = b.thin.normed %>%
  select(-beta, -type_rg) %>%
  spread(Probe, q)

## training model

training = inner_join(b.wide.normed, training.ancestry)

x.training = training %>%
  select(-gsm.id, -ancestry) %>%
  as.matrix

y.training = factor(training$ancestry)

rf = randomForest(x.training, y.training)


## testing model

x.testing = b.wide.normed %>%
  select(-gsm.id) %>%
  as.matrix

p = predict(rf.100, newdata=x.testing)


predicted.ancestry = data.frame(predicted.ancestry=p,
                                gsm.id=b.wide.normed$gsm.id) %>%
  inner_join(sample.info %>%
               select(gsm.id, series.id, ancestry))


## evalution 

acc = predicted.ancestry %>%
  filter(!is.na(ancestry)) %>%
  filter(ancestry %in% c("EUR", "AFR", "ASN")) %>%
  group_by(series.id) %>%
  summarize(acc=sum(ancestry == predicted.ancestry) / n())