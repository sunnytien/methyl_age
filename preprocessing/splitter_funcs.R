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
  
  d = fread(filename)
  
  if(ag("meth")) return(get.beta.meth(d))
  if(ag("signal")) return(get.beta.signal(d))
  stop("Did not detect filetype (methylation or signal)\n")
}

get.beta.meth = function(d){
  
  meth = d %>%
    select(contains("meth", ignore.case=T)) %>%
    select(-contains("unmeth", ignore.case=T)) %>%
    as.matrix %>%
    apply(2, as.numeric)
  
  unmeth = d %>%
    select(contains("unmeth", ignore.case=T)) %>%
    as.matrix %>%
    apply(2, as.numeric)
  
  if(any(dim(meth) != dim(unmeth))) stop("meth and unmeth matrices have different dimensions")
  
  colnames(meth) = gsub(".Meth.+$", "", colnames(meth), ignore.case=T, perl=T)
  colnames(unmeth) = gsub(".Unmeth.+$", "", colnames(meth), ignore.case=T, perl=T)
  
  if(any(colnames(meth) != colnames(unmeth))) stop("Different column names for meth and unmeth")
  
  beta = meth / (meth + unmeth + 100)
  colnames(beta) = colnames(meth)
  
  id.column = grep("ID|Info|V1", colnames(d), ignore.case=T, value=T) %>%
    grep("REF|Targe|IlluminID|geneInfo|V1", ., ignore.case=T, value=T)
  
  ids = d %>%
    select(contains(id.column)) %>%
    as.data.frame
  
  rownames(beta) = ids[[1]]
  return(beta)
  
}

get.beta.signal = function(d){
  signal = d %>%
    select(contains("signal", ignore.case=T))
  
  # methylated is "Signal B"
  meth = signal %>%
    select(ends_with("B", ignore.case=T)) %>%
    as.matrix %>%
    apply(2, as.numeric)
  
  # unmethylated is "Signal A"
  unmeth = signal %>% 
    select(-ends_with("B", ignore.case=T)) %>% 
    as.matrix %>%
    apply(2, as.numeric)
  
  if(any(dim(meth) != dim(unmeth))) stop("meth and unmeth matrices have different dimensions")
  
  colnames(meth) = gsub(".signal.+$", "", colnames(meth), ignore.case=T, perl=T)
  colnames(unmeth) = gsub(".signal.+$", "", colnames(meth), ignore.case=T, perl=T)
  
  if(any(colnames(meth) != colnames(unmeth))) warn("Different column names for meth and unmeth")
  
  beta = meth / (meth + unmeth + 100)
  colnames(beta) = colnames(meth)
  
  id.column = grep("ID|Info|V1$", colnames(d), ignore.case=T, value=T) %>%
    grep("REF|Targe|geneInfo|IlluminID|V1", ., ignore.case=T, value=T)
  
  if(ncol(ids) > 1) warn("Multiple ID columns found\n")
  
  id.column = if(all(id.column == c("V1", "TargetID"))) "TargetID" else
    id.column
  
  ids = d %>%
    select(contains(id.column)) %>%
    as.data.frame
  
  rownames(beta) = ids[[1]]
  return(beta)
}

write.beta = function(gse.id, beta){
  write.probe = function(id){
    
    probe.file = paste(DATA_DIR, gse.id, "/GSM/", id, ".Rdata", sep="")
    if(!file.exists(probe.file)){
      x = beta[,id]
      names(x) = rownames(beta)
    
      save(x, file=probe.file)
    }
    return(T)
  }
  
  sapply(colnames(beta), write.probe)
}

write.mapping = function(gse.id, mapping){
  mapping.file = paste(DATA_DIR, gse.id, "/GSM/mapping.Rdata", sep="")
  save(mapping, file=mapping.file)
}

splitter = function(gse.id){
  cat(paste("Initiating splitting of", gse.id, "\n"))
  
  beta = get.beta(gse.id)
  write.beta(gse.id, beta)
  
  mapping = map.to.gse(gse.id, beta)
  write.mapping(gse.id, mapping)
}

map.to.gse = function(gse.id, beta){
  
  if(file.exists(paste(DATA_DIR, gse.id, "/GSM/mapping.Rdata", sep=""))){
    return(get(load(paste(DATA_DIR, gse.id, "/GSM/mapping.Rdata", sep=""))))
  }
  
  if(gse.id %in% c("GSE56581",
                   "GSE56105",
                   "GSE43414",
                   "GSE53051",
                   "GSE56046",
                   "GSE61431",
                   "GSE58477",
                   "GSE61151",
                   "GSE36064",
                   "GSE41826")) return(map.to.gse.custom(gse.id))
  
  if(all(grepl("GSM[0-9]+", colnames(beta)))) return(map.to.gse.simple(beta))
  
  if(file.exists(paste(DATA_DIR,
                       gse.id, 
                       "/", 
                       gse.id, 
                       "_series_matrix.txt", 
                       sep=""))) return(map.to.gse.cor(gse.id, beta))
  
  return(NA)
}

map.to.gse.simple = function(beta){
  cat("Looking like ids are already GSM ids\n")
  find.id = function(x) str_match(x, "GSM[0-9][0-9][0-9][0-9][0-9][0-9][0-9]")[1]
  data.frame(raw.id=colnames(beta),
             normed.id=sapply(colnames(beta), find.id),
             r=NA,
             stringsAsFactors=F)
}

clean.series.file = function(series.file){
  cleaned.file = paste(series.file, "cleaned", sep="_")
  
  if(!file.exists(cleaned.file)){
    cmd = paste("cat ", series.file, " | grep -v ! >", cleaned.file, sep="")
    system(cmd)
  }
  cleaned.file
}

map.to.gse.cor = function(gse.id, beta){
  
  cat("Mapping via correlation\n")
  
  type1 = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest %>%
    as.data.frame %>%
    mutate(probe.id=row.names(.)) %>%
    filter(Type=="I") %>%
    {.$probe.id}
  
  series.file = list.files(paste(DATA_DIR, gse.id, sep=""), full.names=T) %>%
    grep("series_matrix", ., value=T) %>%
    grep("txt$", ., value=T) %>%
    clean.series.file
  
  d = fread(series.file)
  
  beta.normed = d %>%
    select(starts_with("GSM")) %>%
    as.matrix %>%
    apply(2, as.numeric)
  
  ids = d %>%
    select(contains("ID")) %>%
    as.data.frame
  
  rownames(beta.normed) = ids[[1]]
  
  sites = intersect(rownames(beta.normed), rownames(beta)) %>%
    {
      .[. %in% type1]
    } %>%
    sample(5000)
  
  correlations = cor(beta[sites, ], beta.normed[sites, ], method="spearman", use="pair") %>%
    as.data.frame %>%
    mutate(raw.id=row.names(.)) %>%
    gather(normed.id, r, -raw.id) %>%
    group_by(raw.id) %>%
    filter(r == max(r)) %>%
    ungroup
  
  n.normed = correlations %>%
    select(normed.id) %>%
    distinct %>%
    nrow
  
  n.raw = correlations %>%
    select(raw.id) %>%
    distinct %>%
    nrow
  
  if(n.normed < n.raw) stop("Multiple raw samples mapping to same normalized sample\n")
  if(min(correlations$r) < 0.95) stop("Poor correlation found in samples\n")
  
  return(correlations)
}

map.to.gse.custom = function(gse.id){
  if(gse.id=="GSE56105") return(mapper.GSE56105())
  if(gse.id=="GSE56581") return(mapper.GSE56581())
  if(gse.id=="GSE43414") return(mapper.GSE43414())
  if(gse.id=="GSE53051") return(mapper.GSE53051())
  if(gse.id=="GSE56046") return(mapper.GSE56046())
  if(gse.id=="GSE61431") return(mapper.GSE61431())
  if(gse.id=="GSE58477") return(mapper.GSE58477())
  if(gse.id=="GSE61151") return(mapper.GSE61151())
  if(gse.id=="GSE36064") return(mapper.GSE36064())
  if(gse.id=="GSE41826") return(mapper.GSE41826())
  
  stop("Could not find custom function for", gse.id, "\n")
}

mapper.GSE56105 = function(){
  
  header = paste(DATA_DIR, "GSE56105/GSE56105_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("Sample_title|Sample_geo_accession", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id))
  
  return(correlations)
  
}

mapper.GSE56581 = function(){
  header = paste(DATA_DIR, "GSE56581/GSE56581_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("Sample_title|Sample_geo_accession", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id)) %>%
    mutate(raw.id=str_match(raw.id, "[0-9]+"))
  
}

mapper.GSE43414 = function(){
  header = paste(DATA_DIR, "GSE43414/GSE43414_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("Sample_title|Sample_geo_accession", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id)) %>%
    mutate(Aii=grepl("Aii$", raw.id)) %>%
    mutate(raw.id=str_match(raw.id, "[0-9]+_R[0-9][0-9]C[0-9][0-9]")[,1]) %>%
    mutate(raw.id=ifelse(Aii, paste(raw.id, "1", sep="."), raw.id))
  
  return(correlations)
  
}

mapper.GSE53051 = function(){
  header = paste(DATA_DIR, "GSE53051/GSE53051_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("Sample_description|Sample_geo_accession", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[2]], 
                            normed.id=header[[1]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id))
  
  return(correlations)
  
}

mapper.GSE56046 = function(){
  header = paste(DATA_DIR, "GSE56046/GSE56046_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("Sample_title|Sample_geo_accession", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id)) %>%
    mutate(raw.id=str_match(raw.id, "[0-9]+")[,1])
  
  return(correlations)
  
}

mapper.GSE61431 = function(){
  header = paste(DATA_DIR, "GSE61431/GSE61431_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("barcode|ID_REF", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id)) %>%
    mutate(raw.id=gsub("barcode: ", "", raw.id))
  
  return(correlations)
}

mapper.GSE58477 = function(){
  
  header = paste(DATA_DIR, "GSE58477/GSE58477_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("Sample_title|ID_REF", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id)) %>%
    mutate(raw.id=str_match(raw.id, "[0-9]+_R.+$")[,1])
  
  return(correlations)
  
}

mapper.GSE61151 = function(){
  
  header = paste(DATA_DIR, "GSE61151/GSE61151_series_matrix.txt", sep="") %>%
    readLines %>%
    grep("Sample_title|ID_REF", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id)) %>%
    mutate(raw.id=str_match(raw.id, "p[0-9].+$")[,1])
  
  return(correlations)
  
}

mapper.GSE36064 = function(){
  header = paste(DATA_DIR, "GSE36064/GSE36064_series_matrix.txt", sep="") %>%
    readLines(200) %>%
    grep("Sample_title|ID_REF", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id)) %>%
    mutate(raw.id=str_match(raw.id, "CHB[0-9]+"))
  
  return(correlations)
}

mapper.GSE41826 = function(){
  header = paste(DATA_DIR, "GSE41826/GSE41826_series_matrix.txt", sep="") %>%
    readLines(200) %>%
    grep("Sample_description|ID_REF", ., value=T) %>%
    strsplit("\t")
  
  correlations = data.frame(raw.id=header[[1]], 
                            normed.id=header[[2]],
                            stringsAsFactors=F) %>%
    mutate(raw.id=gsub("\"", "", raw.id)) %>%
    mutate(normed.id=gsub("\"", "", normed.id)) %>%
    filter(grepl("GSM", normed.id))
}