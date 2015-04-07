library("dplyr")
library("data.table")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("tidyr")
library("randomForest")

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

probes = beta %>%
  filter(Probe %in% aimms) %>%
  collect

data = probes %>%
  gather(gsm.id, beta, starts_with("GSM")) %>%
  inner_join(sample.info %>% select(gsm.id, ancestry, tissue)) %>%
  spread(Probe, beta)

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

rf = randomForest(x, y, ntree=2000)
