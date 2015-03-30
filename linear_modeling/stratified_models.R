library("dplyr")
library("lme4")
library("FDb.InfiniumMethylation.hg19")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("tidyr")
library("BatchJobs")
library("magrittr")
library("car")
library("doParallel")

source("./linear_modeling//stratified_model_funcs.R")

load("./data/sample.info.Rdata")
load("./data/predicted.ancestry.Rdata")

select = dplyr::select
group_by = dplyr::group_by
mutate = dplyr::mutate
llply = plyr::llply

#### RUNNING YJ NORMALIZATION ####

yj = powerTransform(sample.info$age, family="yjPower")

sample.info %<>% 
  mutate(age.normed=yjPower(age, coef(yj))) %>%
  mutate(age.normed=scale(age.normed, center=F))

#### FILTERING PROBES ####

## remove X/Y probes
chrXY = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations %>%
  as.data.frame %>%
  mutate(Probe=rownames(.)) %>%
  filter(chr %in% c("chrX", "chrY"))

## remove multiple matching probes
multiple_match = read.table("./data/probes_MM.txt", stringsAsFactors=F) %>%
  select(Probe=V1)

## remove probes with common SNP
snps = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$SNPs.137CommonSingle %>%
  as.data.frame  %>%
  mutate(Probe=rownames(.)) %>%
  filter(!is.na(Probe_maf)) %>%
  filter(Probe_maf > 0.05)

## get remaining probes
db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

probes = beta %>%
  select(Probe) %>%
  collect() %>%
  anti_join(chrXY) %>%
  anti_join(multiple_match) %>%
  anti_join(snps) %>%
  .$Probe

#### MAPPING TO TSS ####

hm450 = get450k()

## map probes to Transcription Start Sites
mapping = getNearestTSS(hm450[probes]) %>%
  as.data.frame %>%
  mutate(Probe=row.names(.)) %>%
  filter(distance<=1500)

## For each gene, choose isoform with the most probes 
## Require at least 4 probe sites
transcripts = mapping %>%
  group_by(nearestTranscript, nearestGeneSymbol) %>%
  summarize(N=n()) %>%
  ungroup %>%
  group_by(nearestGeneSymbol) %>%
  filter(N==max(N)) %>%
  ungroup %>%
  filter(N>=4)

## merge all of the relevant information 
## for all the probes
## need to know which gene, which isoform, and position

probe.info = hm450 %>%
  as.data.frame %>%
  mutate(Probe=rownames(.)) %>%
  select(Probe, Chr=seqnames, Position=probeTarget) %>%
  inner_join(mapping %>% select(Probe, nearestGeneSymbol, nearestTranscript)) %>%
  semi_join(transcripts)

probe.list = probe.info %>%
  split(., .$nearestGeneSymbol)

registerDoParallel(10)

models = llply(probe.list, 
               run.model, 
               sample.info, 
               predicted.ancestry, 
               save=T,
               .parallel=T)
