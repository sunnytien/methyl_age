library("dplyr")
library("nlme")
library("FDb.InfiniumMethylation.hg19")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("magrittr")
library("tidyr")

run.lme = function(probes, gene, beta){
  
  blood = c("Leukocytes",
            "Lymphoblasts",
            "Monocytes",
            "T-cells",
            "Whole Blood")
  
  methylation = beta %>%
    filter(Probe %in% probes) %>%
    collect %>%
    gather(gsm.id, beta, -Probe) %>%
    mutate(M=log(beta / (1-beta)))
  
  data = methylation %>%
    inner_join(sample.info) %>%
    mutate(ancestry=ifelse(is.na(ancestry),
                           sample(c("EUR", "AFR", "ASN"), nrow(methylation), replace=T),
                           ancestry)) %>%
    inner_join(probe.info) %>%
    filter(!is.na(M)) %>%
    filter(!is.na(age)) %>%
    mutate(gsm.id=factor(gsm.id)) %>%
    as.data.frame %>%
    filter(ancestry %in% c("EUR", "AFR", "ASN")) %>%
    filter(tissue != "Lymphocysts") %>%
    mutate(tissue_class=ifelse(tissue %in% blood, "liquid", "solid")) %>%
    mutate(Probe=factor(Probe))
  
  contrasts(data$Probe) = contr.sum(length(levels(data$Probe)))
  
  m.car = lme(M ~ age + tissue_class + ancestry, 
               data=data,
               random=list(~1|gsm.id),
               correlation=corCAR1(0.99, form= ~Position|gsm.id),
               method="REML")
  
  m.lin = lme(M ~ age*Probe, 
               data=data,
               random=list(~1|gsm.id, ~age*Probe|tissue_class),
               method="REML")
  
  
  
  m.exp = lme(M ~ age + tissue_class + ancestry, 
              data=data,
              random=list(~Position|gsm.id),
              correlation=corLin(form= ~Position|gsm.id),
              method="REML")
  
  m.gaus = lme(M ~ age + tissue_class + ancestry + Probe, 
               data=data,
               random=list(~1|gsm.id),
               correlation=corGaus(form= ~Position|gsm.id),
               method="REML")
  
  m = lme(M ~ age*Probe,
          random=list(~1|gsm.id, ~age*Probe|ancestry),
          data=data)
  
}

load("./data/sample.info.Rdata")
load("./data/predicted.ancestry.Rdata")

select = dplyr::select
group_by = dplyr::group_by
mutate = dplyr::mutate

db = src_sqlite("./data/BMIQ.db")
beta = tbl(db, "BMIQ")

probes = beta %>%
  select(Probe) %>%
  collect() %>%
  .$Probe

## remove X/Y probes
## remove probes with common SNP
## remove multiple matching probes

hm450 = get450k()

## map probes to Transcription Start Sites
mapping = getNearestTSS(hm450[probes]) %>%
  as.data.frame %>%
  mutate(Probe=row.names(.)) %>%
  filter(distance<=1500)

## For each gene, choose isoform with the most probes 
## in TSS
transcripts = mapping %>%
  group_by(nearestTranscript, nearestGeneSymbol) %>%
  summarize(N=n()) %>%
  ungroup %>%
  group_by(nearestGeneSymbol) %>%
  filter(N==max(N)) %>%
  ungroup %>%
  filter(N>=5)

## merge all of the relevant information 
## for all the probes
## need to know which gene, which isoform, and position

probe.info = hm450 %>%
  as.data.frame %>%
  mutate(Probe=rownames(.)) %>%
  select(Probe, Chr=seqnames, Position=probeTarget) %>%
  inner_join(mapping %>%
               select(Probe, nearestGeneSymbol, nearestTranscript)) %>%
  semi_join(transcripts)

probes = probe.info %>%
  filter(nearestGeneSymbol=="GRM2") %>%
  .$Probe
