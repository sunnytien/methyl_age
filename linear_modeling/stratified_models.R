library("BatchJobs")

load("./data/horvath_ages.Rdata")
source("./linear_modeling/stratified_model_funcs.R")
source("./linear_modeling/util.R")

probe.infos = get.probe.infos() # definted in util.R

run.model.wrapper = function(probe.info, ...){
  data = get.model.data(probe.info) # get data for given set of probes
  run.model(data, ...) # train lmer
}

# I'm accessing BIC via the BatchJobs interface
# it allows distributed computation in R 
# see https://github.com/tudo-r/BatchJobs for more info

# first we create a registry
# status about all the jobs will be logged into this registry
reg = makeRegistry("models", src.files=c("./linear_modeling/util.R",
                                         "./linear_modeling/stratified_model_funcs.R"))


# then, we apply run.model.wrapper to each element
# of probe.infos
batchMap(reg, 
         run.model.wrapper, 
         probe.infos, 
         more.args=list(horvath_ages=horvath_ages))

# submit jobs to cluster
submitJobs(reg, 
           chunk(findNotSubmitted(reg), 
                 n.chunks=150))