library("wateRmelon")
library("dplyr")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library("BatchJobs")
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)


bmiq.normalization = function(filename, types, subsample=NULL){
  
  require("wateRmelon")
  
  # filename handling
  betas = get(load(filename))
  base = gsub(".Rdata$", "", filename)
  outfile = paste(base, "_BMIQ.Rdata", sep="")
  
  # sanity check
  betas[betas < 0] = 0 
  betas[betas > 1] = 1
  betas = betas[!is.na(betas)]
  betas = betas[!is.null(betas)]
  
  # making sure beta and types match
  betas = betas[names(betas) %in% names(types)] 
  beta.types = types[names(betas)] 
  
  # conduct normalization and save
  betas.normed = BMIQ(betas, beta.types)
  save(betas.normed, file=outfile)

  return(T)
}


## get files for normalization
files = list.files("~/data/methyl_age/GEO/GSM", 
                   recursive=T,
                   full.names=T) %>%
  grep("Rdata$", ., value=T) %>%
  grep("BMIQ", ., value=T, invert=T)

cat(paste("Normalizing", length(files), "files\n"))


## getting types for sites
types = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest$Type %>%
  as.character %>%
  factor %>%
  as.numeric
names(types) = rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Manifest)

## discarding sites in ChrX and ChrY
## BMIQ funtion will automatically discard sites
## for which there is no type
chrXY = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations %>%
  as.data.frame %>%
  mutate(probe=rownames(.)) %>%
  filter(chr %in% c("chrY", "chrX"))

## discarding random probes
random = IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other %>%
  as.data.frame %>%
  mutate(probe=rownames(.)) %>%
  filter(Random_Loci!="")

## discarding multiple matching probes

types = types[!(names(types) %in% chrXY$probe)]
types = types[!(names(types) %in% random$probe)]

## conduct normalization
reg = makeRegistry("BMIQ", packages=c("wateRmelon"))
batchMap(reg, bmiq.normalization, files, more.args=list(types=types))
submitJobs(reg, chunk(findNotSubmitted(reg), n.chunks=50))

