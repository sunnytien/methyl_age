get.header = function(gse.id){
  list.files(paste(DATA_DIR, gse.id, sep=""), full.names=T) %>%
    grep("series_matrix", ., value=T) %>%
    file %>%
    readLines(200) %>%
    grep("^!|ID.REF", ., value=T)
}

gse.mapper = function(gse.id){
  
  gse.id = str_match(filename, "GSE[0-9]+")[1,1]
  
  id = gsub("^.+/GSE[0-9]+_", "", filename) %>%
    gsub(".Rdata$", "", .)
  
  header = get.header(gse.id)
  
  ## first try finding exact match
  f = if(any(grepl(id, header))){ 
    grep
  } else if(any(agrepl(id, header))){
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
  if(length(ind) > 1) stop("Multiple matches for ID found")
  
  data.frame(series.id=gse.id,
             gsm.id=gsms[[ind]],
             raw.id=id)
  
}

