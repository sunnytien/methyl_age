library("BatchJobs")

source("splitter_funcs.R")

gses = list.files(DATA_DIR) %>%
  grep("GSE", ., value=T)

reg = makeRegistry("splitter", 
                   packages=c("dplyr", 
                              "IlluminaHumanMethylation450kanno.ilmn12.hg19",          
                              "tidyr", 
                              "data.table", 
                              "stringr"),
                   src.files="splitter_funcs.R")

batchMap(reg, splitter, gses)
submitJobs(reg)