library("dplyr")
library("nlme")

run.model = function(probe.info, sample.info, predicted.ancestry, db=NULL){
  require("dplyr")
  require("lmer")
  require("tidyr")
  
  liquid_tissues = c("Leukocytes", "Lymphoblasts", "Lymphocytes", "Monocytes",
                     "T-cells", "Whole Blood")
  
  if(is.null(db)) db = src_sqlite("./data/BMIQ.db")
  
  beta = tbl(db, "BMIQ") %>%
    filter(Probe %in% probe.info$Probe) %>%
    collect
  
  beta.thin = beta %>%
    gather(gsm.id, beta, starts_with("GSM")) %>%
    mutate(M=log(beta/(1-beta))) 
  
  data = beta.thin %>%
    inner_join(sample.info) %>%
    inner_join(predicted.ancestry) %>%
    mutate(tissue_state = tissue %in% liquid_tissues) %>%
    mutate(Probe=factor(Probe)) %>%
    filter(tissue != "Lymphoblasts")
  
  m = lmer(M ~ age*Probe + tissue_state + ancestry + (1+age+Probe)|tissue + (1+Probe)|gsm.id,
          data=data,
          REML=F,
          verbose=2)
  
}