run.model = function(probe.info, sample.info, predicted.ancestry, db=NULL){
  require("dplyr")
  require("tidyr")
  require("lme4")
  
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
  
  contrasts(data$Probe) = contr.sum(length(levels(data$Probe)))
  
  lc = lmerControl(check.conv.grad     = .makeCC("stop", tol = 1e-3, relTol = NULL),
                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                    check.conv.hess     = .makeCC(action = "stop", tol = 1e-6),
                   optCtrl=list(maxfun=4e5))
  
  m = lmer(M ~ age*Probe + tissue_state + predicted.ancestry + (1|tissue) + (Probe|tissue) + (1|gsm.id),
          data=data,
          REML=F,
          verbose=2,
          control=lc)
  
  return(m)
}