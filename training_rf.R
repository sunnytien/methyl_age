library("dplyr")
library("tidyr")
library("ggvis")
library("randomForest")

get.ancestry = function(x){
  if(grepl("African", x)) return("AFR")
  if(grepl("Caucasian", x)) return("EUR")
  if(grepl("Han|Japanese",x )) return("ASN")
  return(NA)
}

### Processing data

load("./data/aimms.Rdata")
load("./data/sample.info.Rdata")

training.ancestry = read.delim("./data/training_ancestry.txt",
                               header=T,
                               stringsAsFactors=F) %>%
  mutate(ancestry=sapply(ancestry, get.ancestry)) %>%
  mutate(gsm.id=gsub(" ", "", gsm.id))

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

names(aimms) = NULL

beta.aimms = beta %>%
  filter(Probe %in% aimms) %>%
  collect

b = beta.aimms %>%
  gather(gsm.id, beta, -Probe)

### Training modeling using lymphoblasts

training = b %>%
  inner_join(training.ancestry) %>%
  spread(Probe, beta)

x.training = training %>%
  select(-gsm.id, -ancestry) %>%
  as.matrix

y.training = factor(training$ancestry)

rf = randomForest(x.training, y.training,
                  ntree=2000)

save(rf, file='./data/ancestry_rf.Rdata')

## Applying model to rest of data

x.testing = beta.aimms %>%
  select(-Probe) %>%
  as.matrix %>%
  t
colnames(x.testing) = beta.aimms$Probe

p = predict(rf, newdata=x.testing)

predicted.ancestry = data.frame(predicted.ancestry=p,
                                gsm.id=names(p),
                                stringsAsFactors=F)

sample.info2 = inner_join(sample.info, predicted.ancestry)

tmp = sample.info2 %>%
  filter(!is.na(ancestry))

series = tmp %>% 
  filter(ancestry %in% c("EUR", "AFR", "ASN")) %>%
  group_by(series.id, tissue) %>% 
  summarize(acc=sum(ancestry==predicted.ancestry)/n())

