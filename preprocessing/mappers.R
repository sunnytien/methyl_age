library("BatchJobs")
library("magrittr")

source("mapper_funcs.R")

files = list.files("~/data/methyl_age/GEO/GSM", full.names=T)
# 
# reg = makeRegistry("mapping", 
#                    packages=c("dplyr", "stringr"),
#                    src.files="mapper_funcs.R")
# batchMap(reg, gse.mapper, files)
# 
# submitJobs(reg, 1:10)
# submitJobs(reg, chunk(findNotSubmitted(reg),
#                       n.chunks=30))

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

files = list.files("./data/GSE40279/GSM", full.names=T)
