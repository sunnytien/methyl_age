#DATA_DIR = "~/data/methyl_age/GEO/"
DATA_DIR = "./data/"

get.header = function(gse.id){
  matrix.file = list.files(paste(DATA_DIR, gse.id, sep=""), full.names=T) %>%
    grep("series_matrix.txt$", ., value=T) 
  
  if(length(matrix.file) == 0) stop("No series matrix file found")
  if(length(matrix.file) > 1) stop("Multiple series matrix files found")
  
  matrix.file %>%
    file %>%
    readLines(200) %>%
    grep("^!|ID.REF", ., value=T)
}

gse.mapper = function(filename){
  

  require("dplyr")
  require("stringr")
  
  gse.id = str_match(filename, "GSE[0-9]+")[1,1]
  
  id = gsub("^.+/GSE[0-9]+_", "", filename) %>%
    gsub(".Rdata$", "", .)
  
  header = get.header(gse.id)
  
  ## first try finding exact match
  f = if(any(grepl(id, header))){ 
    grep
  } else if(any(agrepl(id, header))){ # then try fuzzy match
    agrep
  } else{
    stop(paste("Could not find id for: ", gse.id, id, "\n"))
  }
  
  id.line = f(id, header, value=T)
  gsm.line = grep("ID.REF", header, value=T)
  
  ids = id.line %>%
    gsub("\"", "", .) %>%
    gsub("\'", "", .) %>%
    strsplit("\t") %>%
    .[[1]]
  
  gsms = gsm.line %>%
    gsub("\"", "", .) %>%
    gsub("\'", "", .) %>%
    strsplit("\t") %>%
    .[[1]]
  
  ind = f(id, ids)
  
  ## trying to rescue a multiple match here
  ## GSM10 and GSM100 should resolve uniquely
  ## after this statement
  if(length(ind) > 1){
    ind = f(paste(id, "[^0-9]", sep=""), ids)
  }
  
  if(length(ind) > 1) stop(paste("Multiple matches for ID found", gse.id, id))
  if(length(ind) == 0) stop(paste("No matches for ID found", gse.id, id))

  data.frame(series.id=gse.id,
             gsm.id=gsms[[ind]],
             raw.id=id,
             stringsAsFactors=F)
  
}
