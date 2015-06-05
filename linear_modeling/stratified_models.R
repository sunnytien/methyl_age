library("BatchJobs")

load("./data/horvath_ages.Rdata")
source("./linear_modeling/stratified_model_funcs.R")
source("./linear_modeling/util.R")



run.model.wrapper = function(probe.info, ...){
  data = if(is.data.frame(probe.info)){ 
      get.model.data(probe.info)
    } else{
      get(load(probe.info)) 
    }# get data for given set of probes
  run.model(data, ...) # train lmer
}

# I'm accessing BIC via the BatchJobs interface
# it allows distributed computation in R 
# see https://github.com/tudo-r/BatchJobs for more info

# first we create a registry
# status about all the jobs will be logged into this registry
reg = makeRegistry(paste("models", Sys.time(), sep="_") %>% gsub("[ :-]", "_", .), 
                   src.files=c("./linear_modeling/util.R",
                                         "./linear_modeling/stratified_model_funcs.R"))

#probe.infos = get.probe.infos() # definted in util.R 
probe.infos = list.files("./data/model_data", full.names=T) # this directory on bic contains the data.frames already
                                                            # just need to load this files, no need to do SQL calls

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