library("purrr")

get.params = function(method, ...){
  params = method %>%
    purrr::when(
      # regression methods
      . == "glmnet" ~ get.params.glmnet(x, y),
     . == "pcr" ~ get.params.pcr(...),
     . == "pls" ~ get.params.pls(...),
     . == "lda" ~ get.params.lda(...),
      
      # kernel methods
      . == "gaussprLinear" ~ get.params.gaussprLinear(...),
      . == "gaussprRadial" ~ get.params.gaussprRadial(...),
      . == "rvmLinear" ~ get.params.rvmLinear(...),
      . == "rvmRadial" ~ get.params.rvmRadial(...),
      . == "svmRadialCost" ~ get.params.svmRadialCost(...),
      . == "svmLinear" ~ get.params.svmLinear(...),
      
      # tree methods
      . == "rf" ~ get.params.rf(...),
      . == "gbm" ~ get.params.gbm(...),
      . == "cubist" ~ get.params.cubist(...),
        ~ NULL
    )
  params = if(!is.null(params)) as.data.frame(params) else NULL
  return(params)
}

## REGRESSION METHODS
get.params.glmnet = function(x, y){
  
  get.lambda = function(alpha) data.frame(lambda=glmnet(x,y)$lambda)
  
  params = data.frame(alpha=seq(0,1, length.out=11)) %>%
    group_by(alpha) %>%
    do(get.lambda(.$alpha[1])) %>%
    ungroup
  
  return(params)
}
get.params.pcr = function(...) data.frame(ncomp=1:5)
get.params.pls = function(...) data.frame(ncomp=1:5)
get.params.ppr = function(...) data.frame(nterms=1:5)

## KERNAL METHODS
get.params.gaussprLinear = function(...) NULL # no parameters
get.params.gaussprRadial = function(...) NULL # no parameters
get.params.rvmLinear = function(...) NULL # no parameters
get.params.rvmRadial = function(...) NULL # no parameters
get.params.svmRadialCost = function(...) NULL # defaults are ok
get.params.svmLinear = function(...) NULL # defaults are ok

## TREE METHODS
get.params.rf = function(...) NULL # defaults are fine
get.params.cubist = function(...){
  return(expand.grid(neighbors=c(1,5,9), 
                      committees=c(1,100)))
}
get.params.gbm = function(...){
  return(expand.grid(n.trees=c(100,1000,10000,100000),
                    interaction.depth=1:3,
                    shrinkage=0.001))
}