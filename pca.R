library("BatchJobs")
library("pcaMethods")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library("dplyr")

subsetter = function(filename, sites){  
  betas = get(load(filename))$nbeta
  
  id = filename %>%
    (function(x) gsub("^.+GSM/", "", x)) %>%
    (function(x) gsub("_BMIQ.Rdata$", "", x))
  
  missing.sites = sites[!(sites %in% names(betas))]
  missing = rep(NA, length(missing.sites))
  names(missing) = missing.sites
  
  betas = c(betas, missing)
  list(betas=betas[sites],
       id=id)
}

files = list.files("~/data/methyl_age/GEO",
                   full.names=T,
                   recursive=T) %>%
  (function(x) grep("BMIQ", x, value=T)) %>%
  (function(x) grep("Rdata$", x, value=T))

sites = sample(rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest),
               1000)
subset.reg = makeRegistry("subsetting", packages=c("wateRmelon", "dplyr"))

batchMap(subset.reg, subsetter, files, more.args=list(sites=sites))

submitJobs(subset.reg, 1)
submitJobs(subset.reg, chunk(findNotSubmitted(subset.reg),
                      n.chunks=10))

betas = loadResults(subset.reg)
betas = do.call(cbind, betas)

betas = betas[,apply(betas, 2, function(x) all(!is.na(x)))]

p = pca(betas[,1:1000], center=T, scale="none")

scores = p %>%
  scores %>%
  as.data.frame %>%
  mutate(id=rownames(.))

scores %>%
  ggvis(~PC1, ~PC2, key := ~id) %>%
  layer_points() %>%
  add_tooltip(all_values)

