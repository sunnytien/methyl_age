library("wateRmelon")
library("dplyr")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)


bmiq.normalization = function(filename, types, subsample=NULL){
  
  base = gsub(".Rdata$", "", filename)
  
  betas = get(load(filename))
  betas = betas[names(betas) %in% names(types)] # discard unknown sites
  betas[betas < 0] = 0 # sanity check
  betas[betas > 1] = 1
  
  beta.types = types[names(betas)] # get types for this betas
  
  betas.normed = BMIQ(betas, beta.types)
  
  save(betas.normed, file=paste(base, 
                                "_BMIQ.Rdata",
                                sep=""))
  
  if(!is.null(subsample)){
    sites = subsample[subsample %in% names(betas)]
    betas.reduced = betas.normed$nbeta[sites]
    
    save(betas.reduced, 
         file=paste(base,
                    "_BMIQ_reduced.Rdata",
                    sep=""))
  }
  return(T)
}

files = list.files("~/data/methyl_age/GEO", 
                   recursive=T,
                   full.names=T) %>%
  (function(x) grep("GSM", x, value=T)) %>%
  (function(x) grep("Rdata$", x, value=T)) %>%
  (function(x) grep("BMIQ", x, value=T, invert=T))


types = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest$Type %>%
  as.character %>%
  factor %>%
  as.numeric
names(types) = rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest)

subsample = sample(names(types), 1000)

reg = makeRegistry("BMIQ", packages=c("wateRmelon"))
batchMap(reg, bmiq.normalization, files, more.args=list(types=types, 
                                                        subsample=subsample))
submitJobs(reg, 1)
submitJobs(reg, chunk(findNotSubmitted(reg), n.chunks=20))

