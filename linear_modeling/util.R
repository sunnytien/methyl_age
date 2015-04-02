library("dplyr")
library("lme4")
library("FDb.InfiniumMethylation.hg19")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("tidyr")
library("BatchJobs")
library("magrittr")
library("car")
library("doParallel")
library("data.table")


select = dplyr::select
group_by = dplyr::group_by
mutate = dplyr::mutate
llply = plyr::llply

get.probe.infos = function(){

  
  #### FILTERING PROBES ####
  
  ## remove X/Y probes
  chrXY = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations %>%
    as.data.frame %>%
    mutate(Probe=rownames(.)) %>%
    filter(chr %in% c("chrX", "chrY"))
  
  ## remove multiple matching probes
  multiple_match = read.table("./data/probes_MM.txt", stringsAsFactors=F) %>%
    select(Probe=V1)
  
  ## remove probes with common SNP
  snps = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$SNPs.137CommonSingle %>%
    as.data.frame  %>%
    mutate(Probe=rownames(.)) %>%
    filter(!is.na(Probe_maf)) %>%
    filter(Probe_maf > 0.05)
  
  ## get remaining probes
  db = src_sqlite("./data/BMIQ.db")
  beta = tbl(db, "BMIQ")
  
  probes = beta %>%
    select(Probe) %>%
    collect() %>%
    anti_join(chrXY) %>%
    anti_join(multiple_match) %>%
    anti_join(snps) %>%
    .$Probe
  
  #### MAPPING TO TSS ####
  
  hm450 = get450k()
  
  ## map probes to Transcription Start Sites
  mapping = getNearestTSS(hm450[probes]) %>%
    as.data.frame %>%
    mutate(Probe=row.names(.)) %>%
    filter(distance<=1500)
  
  ## For each gene, choose isoform with the most probes 
  ## Require at least 4 probe sites
  transcripts = mapping %>%
    group_by(nearestTranscript, nearestGeneSymbol) %>%
    summarize(N=n()) %>%
    ungroup %>%
    group_by(nearestGeneSymbol) %>%
    filter(N==max(N)) %>%
    ungroup %>%
    filter(N>=4)
  
  ## merge all of the relevant information 
  ## for all the probes
  ## need to know which gene, which isoform, and position
  
  probe.info = hm450 %>%
    as.data.frame %>%
    mutate(Probe=rownames(.)) %>%
    select(Probe, Chr=seqnames, Position=probeTarget) %>%
    inner_join(mapping %>% 
                 select(Probe, nearestGeneSymbol, nearestTranscript, distance)) %>%
    semi_join(transcripts)
  
  probe.list = probe.info %>%
    split(., .$nearestGeneSymbol)
  
  return(probe.list)
}

get.model.data = function(probe.info){
  
  load("./data/sample.info.Rdata")
  load("./data/predicted.ancestry.Rdata")
  
  sample.info %<>% 
    mutate(age.normed=ifelse(age<20,
                             log10((age+1)/21) + 1,
                             (age+1)/21))
  
  liquid_tissues = c("Leukocytes", "Lymphoblasts", "Lymphocytes", "Monocytes",
                     "T-cells", "Whole Blood")
  
  db = src_sqlite(get.db.file("./data/BMIQ.db"))
  
  beta = tbl(db, "BMIQ") %>%
    filter(Probe %in% probe.info$Probe) %>%
    collect
  
  beta.thin = beta %>%
    gather(gsm.id, beta, starts_with("GSM")) %>%
    mutate(M=log(beta/(1-beta))) 
  
  data = beta.thin %>%
    inner_join(probe.info) %>%
    inner_join(sample.info) %>%
    inner_join(predicted.ancestry %>%
                 select(gsm.id, predicted.ancestry)) %>%
    mutate(tissue_state = ifelse(tissue %in% liquid_tissues, "Liquid", "Solid")) %>%
    mutate(Probe=factor(Probe),
           predicted.ancestry=factor(predicted.ancestry),
           tissue_state=factor(tissue_state)) %>%
    filter(tissue != "Lymphoblasts") 
  
  return(data)
}

get.db.file = function(db.file){
  db.tmp = paste(Sys.getenv("TMP"), "/db", sep="")
  if(!file.exists(db.tmp)){ 
    cat(paste("Copying data base to", db.tmp, "\n"))
    file.copy(db.file, db.tmp)
  } else {
    cat(paste("Database found at", db.tmp, "\n"))
  }
  return(db.tmp)
}

write.probe.info = function(){
  probe.infos = get.probe.infos()
  probe.info = rbindlist(probe.infos)
  save(probe.info, file="./data/probe.info.Rdata")
}


