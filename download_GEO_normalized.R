## this script downloads normalized
## GEO sets from GEO
## original idea was to use as a sanity check

download.and.save = function(geo.set){
  geo = getGEO(geo.set, getGPL=F)
  save(geo, file=paste("./data/", geo.set, sep=""))
  return(T)
}

library("GEOquery")
library("dplyr")

load("./data/sample.info.Rdata")

series = sample.info %>%
  select(series.id) %>%
  distinct 

