source("./linear_modeling/stratified_model_funcs.R")
source("./linear_modeling/util.R")

probe.infos = get.probe.infos()

data.reg = makeRegistry("model_data", src.files=c("./linear_modeling/util.R",
                                         "./linear_modeling/stratified_model_funcs.R"))

batchMap(data.reg, get.model.data, probe.infos)
submitJobs(data.reg, 1)
submitJobs(data.reg, 
           chunk(findNotSubmitted(data.reg),n.chunks=100))

model.reg = makeRegistry("models", src.files=c("./linear_modeling/util.R",
                                                "./linear_modeling/stratified_model_funcs.R"))
batchMapResults(data.reg, 
                model.reg, 
                run.model,
                save=T)
submitJobs(model.reg, 
           chunk(findNotSubmitted(model.reg),n.chunks=100))