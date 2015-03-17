library("plyr")
library("dplyr")
library("stringr")

## function for adding data for SQLite database

add.to.sql = function(bmiq.file, db, table.name="beta"){
  
  r = get(load(bmiq.file))
  id = str_match(bmiq.file, "GSM[0-9]+")
  
  df = data.frame(gsm.id=id,
                  probe.id=names(r$nbeta),
                  beta=r$nbeta,
                  stringsAsFactors=F)

  if(table.name %in% db_list_tables(db$con)){
    db_insert_into(db$con, table=table.name, df)
  } else {
    copy_to(db, df, table.name, temporary=F)
  }
  
  return(T)
}
add.to.sql.safe = failwith(NA, add.to.sql)

## script

bmiq.files = list.files("./data",
                        recursive=T,
                        full.names=T) %>%
  (function(x) grep("BMIQ.Rdata$", x, value=T))

db = src_sqlite("./data/BMIQ.db", create=T)

cat(paste("Copying", length(bmiq.files), "files to database\n"))

results = llply(bmiq.files, 
                add.to.sql.safe, 
                db=db,
                table.name="test4",
                .progress="text")
