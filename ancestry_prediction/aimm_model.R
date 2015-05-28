## this script trains a model which predicts contental ancestry from
## Ancestry Informative Methylation Markers (aimms)


library("data.table")
library("randomForest")
library("glmnet")
library("dplyr")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450k.db")

load("./data/sample.info.Rdata")

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

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

## processing AIMMS (Ancestry Informative Methylation Markers)

load("./data/aimms.Rdata")

chrXY = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations %>%
  as.data.frame %>%
  mutate(Probe=rownames(.)) %>%
  filter(chr %in% c("chrX", "chrY"))

multiple_match = read.table("./data/probes_MM.txt") %>%
  select(Probe=V1)

aimms = aimms[!(aimms %in% chrXY$Probe)] # remove probes on X or Y chromosome
aimms = aimms[!(aimms %in% multiple_match$Probe)] # remove multiple matching probes
aimms = aimms[aimms %in% (beta %>% select(Probe) %>% collect %>% .$Probe)] # keep only probes in database

## loading training data for ancestry model

d = fread("./data/GSE36369/GSE36369_series_matrix.txt_cleaned") %>%
  filter(ID_REF %in% aimms) 

# create feature matrix for training model
x = d %>%
  select(-ID_REF) %>%
  as.matrix %>%
  t
colnames(x) = d$ID_REF

# get labels for each sample
# we have to extract this information from the 
# header of the GEO file 
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

m = cv.glmnet(x.cleaned,
              pop.cleaned,
              family="multinomial",
              keep=T,
              nfolds=10)

## apply model to entire dataset

features = colnames(x.cleaned)

testing = beta %>%
  filter(Probe %in% features) %>%
  collect 

x.testing = testing %>%
  select(-Probe) %>%
  as.matrix %>%
  t
colnames(x.testing) = testing$Probe

p = predict(m, 
            s="lambda.1se",
            newx=x.testing,
            type="class")[,1]

predicted.ancestry = data.frame(predicted.ancestry=p,
                                gsm.id=rownames(x.testing),
                                stringsAsFactors=F) %>%
  inner_join(sample.info %>% select(gsm.id, series.id, ancestry))


save(predicted.ancestry, file="./data/predicted.ancestry.Rdata")

## evalution 

# overall accuracy for individuals in ASN, EUR and AFR
acc = predicted.ancestry %>%
  filter(!is.na(ancestry)) %>%
  filter(ancestry %in% c("EUR", "AFR", "ASN")) %>%
  group_by(series.id) %>%
  summarize(acc=sum(ancestry == predicted.ancestry) / n())

# pediatric samples with known ancestry
peds = predicted.ancestry %>%
  filter(series.id=="GSE50759") %>%
  filter(!is.na(ancestry))

# confusion matrix
table(peds$ancestry, peds$predicted.ancestry)
