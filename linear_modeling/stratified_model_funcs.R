run.model = function(data, save=T){
  
  require("dplyr")
  require("tidyr")
  require("lme4")
  require("lmerTest")
  require("magrittr")
  
  
  age.trans = function(x) ifelse(x < 20, 
                                 log((x+1)/21),
                                 (x - 20) / (21))
  
  contrasts(data$Probe) = contr.sum(length(levels(data$Probe)))
  contrasts(data$tissue_state) = contr.sum(length(levels(data$tissue_state)))
  
  data %<>% 
    mutate(age.normed=age.trans(age)) %>%
    filter(series.id != "GSE56105")
  
  lc = lmerControl(check.conv.grad     = .makeCC("stop", tol = 1e-3, relTol = NULL),
                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
                   optCtrl=list(maxfun=4e5))
  cat("Training model\n")
  m = lmer(M ~ age.normed*Probe + age.normed*tissue_state + predicted.ancestry + (1|gsm.id) + (1|tissue),
          data=data,
          control=lc)
  
  cat("Running anova\n")
  a = lmerTest::anova(m, type=3)
  
  cat("Getting Coefficients\n")
  co = coef(summary(m))
  
  if(save){ 
    cat("Saving\n")
    save(co, file=paste("./data/models/", data$nearestGeneSymbol[1], ".Rdata", sep=""))
    save(a, file=paste("./data/anovas/", data$nearestGeneSymbol[1], ".Rdata", sep=""))
    rm(a)
    rm(m)
    rm(data)
    return(T)
  } else return(list(model=m, anova=a))
}