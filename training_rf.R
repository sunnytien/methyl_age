library("dplyr")
library("tidyr")
library("ggvis")
library("randomForest")
library("magrittr")

get.ancestry = function(x){
  if(grepl("African", x)) return("AFR")
  if(grepl("Caucasian", x)) return("EUR")
  if(grepl("Han|Japanese",x )) return("ASN")
  return(NA)
}

### Processing data

load("./data/aimms.Rdata")
load("./data/sample.info.Rdata")
load("./data/rf.Rdata")

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

ancestry_rf = randomForest(x.training, y.training,
                  ntree=2000)

save(rf, file='./data/ancestry_rf.Rdata')

## Applying model to rest of data

x.testing = beta.aimms %>%
  select(-Probe) %>%
  as.matrix %>%
  t
colnames(x.testing) = beta.aimms$Probe

p = predict(old.rf, newdata=x.testing[,rownames(importance(old.rf))],
            type="class")

predicted.ancestry = data.frame(predicted.ancestry=p, stringsAsFactors=F) %>%
  mutate(gsm.id = rownames(.)) 

save(predicted.ancestry, file="./data/predicted.ancestry.Rdata")

predicted.ancestry %>%
  filter(!is.na(ancestry)) %>%
  filter(ancestry %in% c("EUR", "ASN", "AFR")) %>%
  group_by(series.id) %>%
  summarize(acc=sum(ancestry == predicted.ancestry)/n())

## finding appropriate cutoff

pdf("./figures/predicting_ancestry.pdf")
predicted.ancestry %>% 
  filter(!is.na(ancestry)) %>%
  select(gsm.id, ancestry, EUR, AFR, ASN) %>%
  gather(predicted.ancestry, prob, -ancestry, -gsm.id) %>%
  filter(ancestry == predicted.ancestry) %>%
  ggvis(~prob, fill=~ancestry) %>%
  group_by(ancestry) %>%
  layer_histograms() %>%
  add_axis("x", title="Predicted Probability of Correct Ancestry")
dev.off()

ancestry.melted = predicted.ancestry %>% 
  filter(!is.na(ancestry)) %>%
  select(gsm.id, ancestry, EUR, AFR, ASN) %>%
  gather(predicted.ancestry, prob, -ancestry, -gsm.id) 
  

ancestry.melted %>%
  filter(!(ancestry %in% c("AFR", "EUR", "ASN"))) %>%
  group_by(gsm.id) %>%
  filter(prob == max(prob)) %>%
  ggvis(~prob, fill=~ancestry) %>%
  group_by(ancestry) %>%
  layer_histograms()

ancestry.melted %>%
  filter(ancestry == "HSN") %>%
  group_by(gsm.id) %>%
  filter(prob == max(prob)) %>%
  group_by(predicted.ancestry) %>%
  summarize(count=n())
