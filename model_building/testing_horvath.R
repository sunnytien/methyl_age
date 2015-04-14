library("dplyr")
library("tidyr")

load("./data/sample.info.Rdata")
load("./data/predicted.ancestry.Rdata")

age.trans = function(x) ifelse(x < 20, 
                               log((x+1)/21),
                               (x - 20) / (21))

age.antitrans = function(x) ifelse(x < 0,
                                   21*exp(x) - 1,
                                   21*x + 20)

coefficients = read.csv("./data/horvath_coefficients.csv") %>%
  select(Probe=CpGmarker, CoefficientTraining, CoefficientTrainingShrunk)

db = src_sqlite("./data/BMIQ.db")

betas = tbl(db, "BMIQ") %>%
  filter(Probe %in% coefficients$Probe) %>%
  collect

horvath_ages = betas %>%
  gather(gsm.id, beta, -Probe) %>%
  inner_join(coefficients) %>%
  group_by(gsm.id) %>%
  summarize(horvath_age_trans=sum(CoefficientTraining*beta) + 0.69550725) %>%
  mutate(horvath_age=age.antitrans(horvath_age_trans))

performance = horvath_ages %>%
  inner_join(sample.info) %>%
  inner_join(predicted.ancestry %>% select(gsm.id, series.id, predicted.ancestry)) %>%
  mutate(error=age-horvath_age) %>%
  filter(series.id != "GSE56105") %>%
  filter(tissue != "Lymphocytes") %>%
  filter(series.id!="GSE60274") %>%
  group_by(predicted.ancestry) %>%
  summarize(cor(age, horvath_age, use="pair", ))

performance %>%
  filter(abs(error) > 50) %>%
  as.data.frame