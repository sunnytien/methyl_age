run.model = function(data, save=T){
  
  require("dplyr")
  require("tidyr")
  require("lme4")
  require("lmerTest")
  
  contrasts(data$Probe) = contr.sum(length(levels(data$Probe)))
  contrasts(data$tissue_state) = contr.sum(length(levels(data$tissue_state)))
  
  lc = lmerControl(check.conv.grad     = .makeCC("stop", tol = 1e-3, relTol = NULL),
                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
                   optCtrl=list(maxfun=4e5))
  
  m = lmer(M ~ age.normed*Probe + tissue_state + predicted.ancestry + (1|gsm.id) + (age|tissue),
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