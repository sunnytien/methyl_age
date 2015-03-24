library("BatchJobs")

files = list.files("~/data/methyl_age/GEO/GSM", full.names=T)

reg = makeRegistry("mapping", 
                   packages=c("dplyr", "stringr"),
                   src.files="mapper_funcs.R")
batchMap(reg, gse.mapper, files)
