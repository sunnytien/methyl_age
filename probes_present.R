## script to count the number of probes present in BMIQ files

library("dplyr")

load("./data/sample.info.Rdata")

get.number.of.probes = function(gsm){
  filename = grep(gsm, gsm.files, value=T)
  if(length(filename) == 0) return(NA)
  
  b = get(load(filename))
  length(b$nbeta)
}

gsm.files = list.files("./data/",
                       recursive=T,
                       full.names=T) %>%
  (function(x) grep("_BMIQ.Rdata$", x, value=T)) 


probes.present = sample.info %>%
  group_by(gsm.id) %>%
  do(get.number.of.probes(.$gsm.id[1]))

save(probes.present, file="./data/probes.present.Rdata")