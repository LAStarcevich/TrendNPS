logit <- function(x){
  y <- log(x / (1 - x))
  return(y)
}

expit <- function(x) {
  y <- 1 / (1 + exp(-x))
  return(y)
}