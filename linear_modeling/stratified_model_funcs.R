## this function runs a linear mixed effect region on a promoter

run.model = function(data, horvath_ages, save=T){
  
  require("dplyr")
  require("tidyr")
  require("lme4")
  require("lmerTest")
  require("magrittr")
  
  m2beta = function(m) exp(m) / (1 + exp(m))
  
  data %<>% mutate(ancestry=factor(ancestry))
  
  # set contrast coding for all catagorical variables
  contrasts(data$Probe) = contr.sum(length(levels(data$Probe)))
  contrasts(data$tissue_state) = contr.sum(length(levels(data$tissue_state)))
  contrasts(data$ancestry) = contr.sum(length(levels(data$ancestry)))
  
  # removing bad data 
  data2 = data %>% 
    filter(series.id != "GSE56105") %>% # this GEO series is biased
    inner_join(horvath_ages) %>%  
    filter(horvath_age < 120)  # if horvath_age > 120, sample is likely cancereous
                              # a couple cancer samples may be lurking in the datasets
  
  # set control so that job crashes if lmer doesn't converge
  # default just throws a warning -- this is more strict
  lc = lmerControl(check.conv.grad     = .makeCC("stop", tol = 1e-3, relTol = NULL),
                    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
                    check.conv.hess     = .makeCC(action = "warning", tol = 1e-6),
                   optCtrl=list(maxfun=4e5))
  
  cat("Training model\n")
  m = lmer(M ~ age.normed*Probe + age.normed*tissue_state + predicted.ancestry + (1|gsm.id) + (1|tissue),
          data=data2,
          control=lc)
  
  cat("Running anova\n")
  a = lmerTest::anova(m, type=3)
  
  cat("Getting Coefficients\n")
  co = m %>%
    summary %>%
    coef %>%
    as.data.frame(stringsAsFactors=F) %>%
    mutate(variable=row.names(.))

  age = co$Estimate[co$variable == "age"]
  intercept = co$Estimate[co$variable=="(Intercept)"]
  tissue_state = co$Estimate[co$variable=="tissue_state1"]
  
  asn = co$Estimate[co$variable=="predicted.ancestryASN"]
  eur = co$Estimate[co$variable=="predicted.ancestryEUR"]
  
  # here, I'm calculating the difference in beta-values caused by aging, ancestry, etc
  # this is a little tricky because we're modeling M-values
  # and there is a non-linear relationship between m-values and beta-values
  # see m2beta for the transformation
  
  # difference in beta-value between 80 year old and 20 year old
  dbeta.age = m2beta(age*80 + intercept) - m2beta(age*20 + intercept)
  
  # difference in beta-value between solid and liquid tissues for a 20 year-old
  dbeta.tissue = m2beta(age*20 + intercept - tissue_state) - m2beta(age*20 + intercept + tissue_state)
  
  # difference between ASN and AFR
  dbeta.asn =  m2beta(age*20 + intercept + asn) - m2beta(age*20 + intercept)
  # difference between EUR and AFR
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