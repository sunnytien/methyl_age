library("dplyr")
library("data.table")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("tidyr")
library("randomForest")
library("magrittr")
library("glmnet")

fst = fread("./data/FstGLOB_CEU_u_YRI_u_CHB.whole_genome.pvalues") %>%
  filter(pvalue < 0.05)

snps = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$SNPs.137CommonSingle %>%
  as.data.frame %>%
  mutate(Probe=row.names(.)) %>%
  mutate(snpID=CpG_rs) %>%
  as.data.table

fst.snps = inner_join(fst, snps)

## trial run

load("./data/sample.info.Rdata")

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

names(aimms) = NULL

probes = beta %>%
  filter(Probe %in% aimms) %>%
  collect

data = probes %>%
  gather(gsm.id, beta, starts_with("GSM")) %>%
  spread(Probe, beta) %>% 
  inner_join(sample.info %>% select(gsm.id, ancestry, tissue))

## training rf

training = data %>%
  filter(tissue == "Lymphoblasts") %>%
  select(-ancestry) %>%
  inner_join(pop.mapping)

x = training %>%
  select(starts_with("cg")) %>%
  as.matrix

y = training$ancestry %>%
  factor

m = cv.glmnet(x, y, 
              family="multinomial",
              alpha=0.5)

## testing rf

testing = data %>%
  filter(!is.na(ancestry)) %>%
  filter(ancestry %in% c("ASN", "AFR", "EUR"))

x = testing %>%
  select(starts_with("cg")) %>%
  as.matrix

y = testing$ancestry

p = predict(m, newx=x, s="lambda.min") 

table(y, p[,1])
plot(m)
help(predict.cv.glmnet)
