source("./linear_modeling/stratified_model_funcs.R")
source("./linear_modeling/util.R")

probe.infos = get.probe.infos()

registerDoParallel(10)

models = llply(probe.infos, 
               run.model,  
               save=T,
               .parallel=T)
