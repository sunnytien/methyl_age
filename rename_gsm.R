## this converts from each set's internal ID
## to GSM ids
## it renames all the files in the GSM directories

rename.gse = function(gse.id){
  
  rename.gsm = function(gsm){
    id = gsm %>%
      gsub(".Rdata$", "", .) %>%
      gsub("^.+/", "", .)
    
    id2 = if(grepl("_meth$", id)) gsub("_meth", "", id) else
      id
    
    #cat(paste("Processing", id, "\n"))
    if(!grepl("GSM[0-9][0-9][0-9][0-9][0-9][0-9]+$", id)){
      if(id2 %in% names(map)){id 
        base = gsub(paste(id, ".Rdata$", sep=""), "", gsm)
        new.file = paste(base, map[id2], ".Rdata", sep="")
        file.rename(gsm, new.file)
      } else{
        warning(paste("COULD NOT FIND MAP FOR", id))
      }
    }
  }
  
  mapping.file = paste(DATA_DIR, gse.id, "/GSM/mapping.Rdata", sep="")
  if(file.exists(mapping.file)){
    cat(paste("Found mapping at", mapping.file,"\n"))
    mapping = get(load(mapping.file))
    
    gsms = list.files(paste(DATA_DIR, gse.id, "/GSM/", sep=""),
                            full.names=T) %>%
      grep("Rdata$", ., value=T) %>%
      grep("mapping", ., value=T, invert=T)
    
    map = mapping$normed.id
    names(map) = mapping$raw.id
    
    cat(paste("Found", length(gsms), "GSM files\n"))
    sapply(gsms, rename.gsm)
    
  } else{
    stop("Could not find mapping file")
  }
}


DATA_DIR = "~/data/methyl_age/GEO/"

gses = list.files(DATA_DIR) %>%
  grep("GSE", ., value=T)

lapply(gses, rename.gse)