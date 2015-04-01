source("./linear_modeling/stratified_model_funcs.R")
source("./linear_modeling/util.R")

probe.infos = get.probe.infos()

reg = makeRegistry("models", src.files=c("./linear_modeling/util.R",
                                         "./linear_modeling/stratified_model_funcs.R"))

batchMap(reg, run.model, probe.infos)
submitJobs(reg, 1)

submitJobs(reg, chunk(findNotSubmitted(reg),
                      n.chunks=150))