## rewriting splitting function to be generalized
## idea is for it to be broadly applicable do all GSE sets

## this turned into a mess
## but it works

library("dplyr")
library('stringr')
library("GEOquery")
library("data.table")
library("tidyr")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("BatchJobs")

#DATA_DIR = "./data/"
DATA_DIR = "~/data/methyl_age/GEO/"

data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

get.beta = function(gse.id){
  cat(paste("Getting beta-values for", gse.id, "\n"))
  
  if(file.exists(paste(gse.id, "/", "idat", sep=""))) get.beta.idat(gse.id) else
    get.beta.notidat(gse.id)
}

get.beta.idat = function(gse.id){
  
  filename = paste(DATA_DIR, gse.id, "/idat", sep="")
  cat(paste("Found idat folder at:", filename, "\n"))

  rg = read.450k.exp(filename, verbose=T)
  meth = preprocessIllumina(rg, bg.correct=T, normalize="no")
  
  beta = getBeta(meth, offset=100)
  
}


get.beta.notidat = function(gse.id){
  
  ag = function(x) any(grepl(x, colnames(d), ignore.case=T))
  
  filename = gse.id  %>%
    paste(DATA_DIR, ., sep="") %>%
    list.files(full.names=T) %>%
    grep(paste(gse.id, gse.id, sep="/"), ., value=T) %>%
    grep("series_matrix", ., value=T, invert=T)
  
  if(length(filename) != 1) stop("Could not locate datafile")
  
  cat(paste("Found possible datafile at:", filename, "\n"))
  
  filename %>%
    fread %>%
    get.beta
}

get.data.file = function(gse.id){
  
  filename = gse.id  %>%
    paste(DATA_DIR, ., sep="") %>%
    list.files(full.names=T) %>%
    grep(paste(gse.id, gse.id, sep="/"), ., value=T) %>%
    grep("series_matrix", ., value=T, invert=T)
  
  if(length(filename) != 1) stop("Could not locate datafile")
  
  cat(paste("Found possible datafile at:", filename, "\n"))
  
  return(filename)
}

get.beta.notidat = function(gse.id){
  
  filename = get.data.file(gse.id)
  
  d = fread(filename)
  
  meth = get.meth(d)
  unmeth = get.unmeth(d)
  if(any(colnames(meth) != colnames(unmeth))) warning("Column names aren't matching")
  
  beta = (meth + 50) / (meth + unmeth + 100)
  
  return(beta)
  
}

get.meth = function(d){
  cat("Finding methylated ")
  meth = if(any(grepl("meth", colnames(d), ignore.case=T))){
    d %>%
      select(contains("meth", ignore.case=T)) %>%
      select(-contains("unmeth", ignore.case=T)) %>%
      as.matrix %>%
      apply(2, as.numeric)
  } else{
    d %>%
      select(contains("signal")) %>%
      select(ends_with("B", ignore.case=T)) %>%
      as.matrix %>%
      apply(2, as.numeric)
  }
  
  rownames(meth) = get.probe.ids(d)
  colnames(meth) = get.colnames(meth)
  return(meth)
}

get.unmeth = function(d){
  unmeth = if(any(grepl("meth", colnames(d), ignore.case=T))){
    d %>%
      select(contains("unmeth", ignore.case=T)) %>%
      as.matrix %>%
      apply(2, as.numeric)
  } else{
    d %>%
      select(contains("signal")) %>%
      select(-ends_with("B", ignore.case=T)) %>%
      as.matrix %>%
      apply(2, as.numeric)
  }
  
  rownames(unmeth) = get.probe.ids(d)
  colnames(unmeth) = get.colnames(unmeth)
  
  return(unmeth)
}

get.probe.ids = function(d){
  tmp = d %>%
    as.data.frame
  
  ## look in the first ten columns for something with a lot
  ## of strings that start with "cg"
  id.column = sapply(tmp[1:10], function(x) sum(grepl("cg", x))) %>%
    which.max
  
  return(tmp[[id.column]])
}

get.colnames = function(m){
  m %>% 
    colnames %>%
    gsub("signal.*$", "", ., ignore.case=T) %>%
    gsub("unmeth.*$", "", ., ignore.case=T) %>%
    gsub("meth.*$", "", ., ignore.case=T) %>%
    gsub(".$", "", .)
}

write.beta = function(gse.id, beta){
  
  write.gsm = function(id){
    x = beta[,id]
    save(x, file=paste(DATA_DIR, "GSM/", gse.id, "_", id, ".Rdata", sep=""))
    return(T)
  }
  
  sapply(colnames(beta), write.gsm)
}

splitter = function(gse.id){
  cat(paste("Initiating splitting of", gse.id, "\n"))
  
  beta = get.beta(gse.id)
  write.beta(gse.id, beta)
}




