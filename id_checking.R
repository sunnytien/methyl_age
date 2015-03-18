library("data.table")
library("plyr")
library("dplyr")
library("tidyr")
library("doMC")

load("./data/sample.info.Rdata")


find.matches = function(gse.id){

  read.series.matrix = function(x){
    
    tmp = tempfile()
    system(paste("cat ", x, " | grep -v ! >", tmp, sep=""))
    
    fread(tmp, na.strings=c("NA", "na", "'null'", "NULL"))
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
  
  
  cat(paste("\nInitiating check for", gse.id, "\n"))
  
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
    as.matrix %>%
    apply(2, as.numeric)
  rownames(beta.normed) = d$ID_REF
  
  beta.raw = meth / (meth + unmeth + 100)
  
  sites = intersect(rownames(beta.normed), rownames(beta.raw)) %>%
    sample(5000)
  
  beta.normed = beta.normed[sites, ]
  beta.raw = beta.raw[sites, ]
  
  
  r = cor(beta.normed, beta.raw, use="pair", method="pearson") %>%
    as.data.frame %>%
    mutate(normed.id=rownames(.)) %>%
    gather(raw.id, r, -normed.id) %>%
    group_by(normed.id) %>%
    filter(r == max(r)) %>%
    ungroup  
  
  if(all(r$normed.id == r$raw.id)) paste(cat(gse.id, "ids match!\n")) else
    cat(paste(gse.id, "ids DO NOT MATCH\n"))
}

find.matches.safe = failwith(NA, find.matches)

series = sample.info %>%
  select(series.id) %>%
  distinct

empty = lapply(series$series.id[2], find.matches.safe)


