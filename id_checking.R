library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("doMC")

load("./data/sample.info.Rdata")

read.series.matrix = function(x){
  
  tmp = tempfile()
  system(paste("cat ", x, " | grep -v ! >", tmp, sep=""))

  fread(tmp)
}

get.file = function(x){
  files = list.files("./data", recursive=T, full.names=T) %>%
                  grep(x, ., value=T)
  function(y) grep(y, files, value=T)
}

get.mm.file = get.file("meth_matrices.Rdata")
get.sm.file = get.file("series_matrix.txt$")

download.sm.file = function(x){
  system("wget")
}

find.matches = function(gse.id){

  mm.file = get.mm.file(gse.id)
  sm.file = get.sm.file(gse.id)
  
  
  if(length(mm.file) == 0) stop(paste("No meth_matrix file found for", gse.id, "\n"))
  if(length(sm.file) == 0) stop(paste("No series matrix file found for", gse.id, "\n"))
  
  cat(paste("Meth matrix file found at:", mm.file, "\n"))
  cat(paste("Series matrix file found at:", sm.file, "\n"))
  
  load("./data/sample.info.Rdata")
  load(mm.file)
  d = read.series.matrix(sm.file)
  
  beta.normed = d %>%
    select(starts_with("GS")) %>%
    as.matrix
  rownames(beta.normed) = d$ID_REF
  
  beta.raw = meth / (unmeth + meth + 100)
  
  sites = intersect(rownames(beta.normed), rownames(beta.raw))
  if(length(sites) < 5e5) stop(paste("Insufficient probes found for", gse.id, "\n")) 
  
  beta.normed = beta.normed[sites, ]
  beta.raw = beta.raw[sites, ]
  
  colnames(beta.normed) = paste(colnames(beta.normed),
                                "normed",
                                sep="_")
  
  colnames(beta.raw) = paste(colnames(beta.raw),
                             "raw",
                             sep="_")
  
  cor(beta.normed, beta.raw) %>%
    as.data.frame %>%
    mutate(id1=rownames(.)) %>%
    gather(id2, r, -id1) %>%
    group_by(id2) %>%
    filter(r == max(r)) %>%
    ungroup  
}

find.matches.safe = failwith(NA, find.matches)

series = sample.info %>%
  select(series.id) %>%
  distinct

registerDoMC(2)

matches = llply(series$series.id, find.matches.safe, .parallel=T)
save(matches, file="./data/matches.Rdata")


