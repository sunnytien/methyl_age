library("BatchJobs")

load("./data/horvath_ages.Rdata")

source("./linear_modeling/stratified_model_funcs.R")
source("./linear_modeling/util.R")

write.model.data = function(probe.info){
  data = get.model.data(probe.info)
  save(data, 
       file=paste("./data/model_data/", probe.info$nearestGeneSymbol[1], ".Rdata", sep=""))
}

run.models2 = function(data.file){
  load(data.file)
  run.model(data)
}

model.data.files = list.files("./data/model_data", full.names=T) %>%
  grep(".Rdata$", ., value=T)

reg = makeRegistry("models2", src.files=c("./linear_modeling/util.R",
                                          "./linear_modeling/stratified_model_funcs.R"))

batchMap(reg, run.models2, model.data.files, more.args=list(horvath_ages=horvath_ages))

submitJobs(reg, chunk(findNotSubmitted(reg), n.chunks=150))