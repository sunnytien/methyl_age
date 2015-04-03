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