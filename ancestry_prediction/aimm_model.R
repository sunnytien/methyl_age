library("data.table")
library("randomForest")
library("glmnet")
library("dplyr")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450k.db")

population = function(x){
  if(grepl("caucasian", x, ignore.case=T)) return("EUR")
  if(grepl("african", x, ignore.case=T)) return("AFR")
  if(grepl("chinese", x, ignore.case=T)) return("ASN")
  if(grepl("japanese", x, ignore.case=T)) return("ASN")
  if(grepl("chimp", x, ignore.case=T)) return("CHIMP")
  NA
}

top_n_sites = function(rf, n){
  importance(rf) %>%
    as.data.frame %>%
    mutate(Probe=row.names(.)) %>%
    arrange(-MeanDecreaseGini) %>%
    head(n) %>%
    .$Probe
}

## preparing AIMMS

load("./data/aimms.Rdata")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

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
aimms = aimms[aimms %in% (beta %>% select(Probe) %>% collect %>% .$Probe)]

## loading data

d = fread("./data/GSE36369/GSE36369_series_matrix.txt_cleaned") %>%
  filter(ID_REF %in% aimms) 

unique.probes = d %>%
  mutate(ID_REF=factor(ID_REF)) %>%
  group_by(ID_REF) %>%
  summarize(N=n()) %>%
  ungroup %>%
  filter(N==1)

d2 = d %>%
  semi_join(unique.probes)

x = d2 %>%
  select(-ID_REF) %>%
  as.matrix %>%
  t
colnames(x) = d2$ID_REF

pop = readLines(file("./data/GSE36369/GSE36369_series_matrix.txt"), 200) %>%
  (function(x) grep("Sample_title", x, value=T, ignore.case=T)) %>%
  (function(x) strsplit(x, "\t")[[1]]) %>%
  (function(x) x[2:length(x)]) %>%
  sapply(population)

names(pop) = NULL

## filtering for bad markers, chimpanzee

rows = pop != "CHIMP"
x.cleaned = x[rows, ]
pop.cleaned = factor(pop[rows])

cols = apply(x.cleaned, 2, function(x) all(!is.na(x)))
x.cleaned  = x.cleaned[,cols]

## training model

set.seed(123)

rf = randomForest(x.cleaned,
                  factor(pop.cleaned),
                  ntree=2000)

sites.1000 = top_n_sites(rf, 1000)

rf.1000 = randomForest(x.cleaned[,sites.1000],
                       factor(pop.cleaned),
                       ntree=2000)

sites.100 = top_n_sites(rf.1000, 100)

rf.100 = randomForest(x.cleaned[,sites.100],
                       factor(pop.cleaned),
                       ntree=2000)


save(rf, rf.1000, rf.100, file="./data/rf.Rdata")


## performance testing

features = rownames(importance(rf))

testing = beta %>%
  filter(Probe %in% features) %>%
  collect 

x.testing = testing %>%
  select(-Probe) %>%
  as.matrix %>%
  t
colnames(x.testing) = testing$Probe

p = predict(rf.100, newdata=x.testing)

predicted.ancestry = data.frame(predicted.ancestry=p,
                                gsm.id=rownames(x.testing),
                                stringsAsFactors=F) %>%
  inner_join(sample.info %>% select(gsm.id, series.id, ancestry))

save(predicted.ancestry, file="./data/predicted.ancestry.Rdata")

## evalution 

acc = predicted.ancestry %>%
  filter(!is.na(ancestry)) %>%
  filter(ancestry %in% c("EUR", "AFR", "ASN")) %>%
  group_by(series.id) %>%
  summarize(acc=sum(ancestry == predicted.ancestry) / n())
