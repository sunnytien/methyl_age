library("dplyr")

train.wrapper = function(method, x, y, ...){
  train(x, y, method=method, tuneGrid=get.params(method, x, y), ...)
}

train.wrapper.safe = failwith(NA, train.wrapper)