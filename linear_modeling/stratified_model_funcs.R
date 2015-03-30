run.model = function(probe.info, sample.info, predicted.ancestry, db=NULL, save=F){
  
  require("dplyr")
  require("tidyr")
  require("lme4")
  require("lmerTest")
  
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
    mutate(tissue_state = ifelse(tissue %in% liquid_tissues, "Liquid", "Solid")) %>%
    mutate(Probe=factor(Probe),
           predicted.ancestry=factor(predicted.ancestry),
           tissue_state=factor(tissue_state)) %>%
    filter(tissue != "Lymphoblasts") 
  
  contrasts(data$Probe) = contr.sum(length(levels(data$Probe)))
  contrasts(data$tissue_state) = contr.sum(length(levels(data$tissue_state)))
  contrasts(data$predicted.ancestry) = contr.sum(length(levels(data$predicted.ancestry)))
  
  lc = lmerControl(check.conv.grad     = .makeCC("stop", tol = 1e-3, relTol = NULL),
                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
                   optCtrl=list(maxfun=4e5))
  
  m = lmer(M ~ age.normed*Probe + age.normed*tissue_state + age.normed*predicted.ancestry + (1 | tissue) + (0 + age|tissue) + (1|gsm.id),
          data=data,
          control=lc)
  
  a = lmerTest::anova(m, type=3)
  
  if(save){ 
    save(m, file=paste("./data/models/", probe.info$nearestGeneSymbol[1], ".Rdata", sep=""))
    save(a, file=paste("./data/anovas/", probe.info$nearestGeneSymbol[1], ".Rdata", sep=""))
    rm(a)
    rm(m)
    rm(data)
    rm(beta.thing)
    rm(beta)
    return(T)
  } else return(list(model=m, anova=a))
}