run.model = function(data, horvath_ages, save=T){
  
  require("dplyr")
  require("tidyr")
  require("lme4")
  require("lmerTest")
  require("magrittr")
  
  m2beta = function(m) exp(m) / (1 + exp(m))
  
  
  data %<>% mutate(ancestry=factor(ancestry))
  
  contrasts(data$Probe) = contr.sum(length(levels(data$Probe)))
  contrasts(data$tissue_state) = contr.sum(length(levels(data$tissue_state)))
  contrasts(data$ancestry) = contr.sum(length(levels(data$ancestry)))
  
  data2 = data %>% 
    filter(series.id != "GSE56105") %>%
    filter(age > 20) %>%
    inner_join(horvath_ages) %>%  
    filter(horvath_age < 100) %>%
    filter(abs(horvath_age - age) < 20)
  
  lc = lmerControl(check.conv.grad     = .makeCC("stop", tol = 1e-3, relTol = NULL),
                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
                   optCtrl=list(maxfun=4e5))
  
  cat("Training model\n")
  m = lmer(M ~ age*Probe + age*tissue_state + age*predicted.ancestry + (1|gsm.id) + (1|tissue),
          data=data2,
          control=lc)
  
  cat("Running anova\n")
  a = lmerTest::anova(m, type=3)
  
  cat("Getting Coefficients\n")
  co = m %>%
    summary %>%
    coef %>%
    as.data.frame %>%
    mutate(variable=row.names(.))

  age = co$Estimate[co$variable == "age"]
  intercept = co$Estimate[co$variable=="(Intercept)"]
  tissue_state = co$Estimate[co$variable=="tissue_state1"]
  
  asn = co$Estimate[co$variable=="predicted.ancestryASN"]
  eur = co$Estimate[co$variable=="predicted.ancestryEUR"]
  
  dbeta.age = m2beta(age*80 + intercept) - m2beta(age*20 + intercept)
  dbeta.tissue = m2beta(age*20 + intercept - tissue_state) - m2beta(age*20 + intercept + tissue_state)
  
  dbeta.asn =  m2beta(age*20 + intercept + asn) - m2beta(age*20 + intercept)
  dbeta.eur =  m2beta(age*20 + intercept + eur) - m2beta(age*20 + intercept)
  
  dbeta = data.frame(variable=c("age", "tissue_state1", "predicted.ancestryASN", "predicted.ancestryEUR"),
                     dBeta=c(dbeta.age, dbeta.tissue, dbeta.asn, dbeta.eur))
  
  co %<>% full_join(dbeta)
  
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