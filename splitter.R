library("BatchJobs")
library("minfi")
library("dplyr")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

process.mset = function(filename){
  
  base = paste(gsub("meth_matrices.Rdata$", "", filename), "GSM/", sep="")
  if(!file.exists(base)) system(paste("mkdir", base))
  
  get.path = function(id) paste(base, id, ".Rdata", sep="")
  
  load(filename)
  
  betas = meth / (meth + unmeth + 100) %>%
    as.data.frame
  
  lapply(names(betas), function(id){
    b = betas[[id]]
    names(b) = rownames(betas)
    save(b, file=get.path(id))
  })
  
  return(rownames(meth))
}

files = list.files("~/data/methyl_age/GEO",
                   full.names=T,
                   recursive=T) %>%
  (function(x) grep("/meth_matrices.Rdata$", x, value=T))

reg = makeRegistry("splitting", packages=c("dplyr", "minfi"))
batchMap(reg, process.mset, files)
submitJobs(reg)

waitForJobs(reg)

annotated.sites = rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest)

sites = loadResults(reg)

unannotated.sites = lapply(sites, function(x) x[!(x %in% annotated.sites)]) %>%
  unlist %>%
  unique