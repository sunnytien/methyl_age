library("dplyr")
library("purrr")
library("stringr")
library("RSQLite")

print(sessionInfo())

addNA = function(x, names){
  
  missing.names = names[!(names %in% names(x))]
  if(length(missing.names) == 0){
    return(x[sort(names(x))])
  } else {
    y = c(x, rep(NA, length(missing.names)))
    names(y) = c(names(x), missing.names)
    y = y[sort(names(y))]
    return(y)
  }
}


db.file = "~/data/methyl_age/GEO/BMIQ.db"

files = list.files("~/data/methyl_age/GEO",
                   full.names=T,
                   recursive=T) %>%
  grep("Rdata", ., value=T) %>%
  grep("BMIQ", ., value=T) %>%
  grep("GSM", ., value=T) 

#ßfiles = sample(files, 1000)

betas = lapply(files, function(x) get(load(x))$nbeta)

names(betas) = str_match(files, "GSM[0-9]+")

sites = lapply(betas, names) %>%
  reduce(union)

cat(paste(length(sites), "sites found\n"))
cat("Adding missing sites to beta vectors\n")

betas2 = lapply(betas, addNA, sites)

rm(betas)

cat("Binding into giant data.frame\n")

betas.df = do.call(data.frame, betas2) %>%
  mutate(Probe=names(betas2[[1]])) %>%
  select(Probe, starts_with("GSM")) %>%
  as.data.frame

rm(betas2)

unlink(db.file)
db = src_sqlite(db.file,
                create=T)

cat("Copying to SQL\n")
betas.sql = copy_to(db, 
                    betas.df,
                    name="BMIQ",
                    temporary=F,
                    indexes=list("Probe"))
