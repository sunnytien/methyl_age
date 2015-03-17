## this script downloads normalized
## GEO sets from GEO
## original idea was to use as a sanity check

library("GEOquery")
library("plyr")
library("dplyr")
library("doMC")

download.and.save = function(geo.set){
  geo = getGEO(geo.set, getGPL=F)
  
  d = dataTable(geo)@table
  save(d, file=paste("./data/GSM/",
                       geo.set, 
                       ".Rdata", 
                       sep=""))
  return(T)
}
download.and.save.safe = failwith(F, download.and.save)

load("./data/sample.info.Rdata")

gsm = sample.info %>%
  select(gsm.id) %>%
  distinct 

registerDoMC(5)

llply(gsm$gsm.id,
      download.and.save.safe,
      .parallel=T)
