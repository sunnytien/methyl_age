library("BatchJobs")
library("magrittr")

source("mapper_funcs.R")

files = list.files("~/data/methyl_age/GEO/GSM", full.names=T)

gse.mapper.safe = failwith(data.frame(series.id=NA,
                                      gsm.id=NA,
                                      raw.id=NA,
                                      stringsAsFactors=F), gse.mapper)

mapping = files %>%
  data.frame(filename=., 
             stringsAsFactors=T) %>%
  group_by(filename) %>%
  do(gse.mapper.safe(.$filename)) %>%
  filter(!is.na(series.id))