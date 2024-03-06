
mydbinom <- function(x, size, prob, log = FALSE){
  if(log){
    ret <- sum(dbinom(x, size = size, prob = prob, log = TRUE))
  } else {

    ret <- prod(dbinom(x, size = size, prob = prob, log = FALSE))
  }
  return(ret)
}

myrbinom <- function(n, size, prob){
  ret <- rbinom(n, size = size, prob = prob)
  return(ret)
}

dmybinom <- nimble::nimbleRcall(
  prototype = function(
    x=double(1), #x is a matrix
    size = double(0),
    prob = double(0),
    log = integer(0, default = 1)
  ) {},
  returnType = double(0), # outcome is a vector
  Rfun = 'mydbinom'
)

rmybinom <- nimble::nimbleRcall(
  prototype = function(
    n = double(0),
    size = double(0),
    prob = double(0)
  ) {},
  returnType = integer(1),
  Rfun = 'myrbinom'
)


mydpois <- function(x, lambda, log = FALSE){
  if(log){
    ret <- sum(dpois(x, lambda, log = TRUE))
  } else {
    ret <- prod(dpois(x, lambda, log = FALSE))
  }
  return(ret)
}

myrpois<- function(n, lambda, lowerBound, includeLowerBound = FALSE){
  if(includeLowerBound){
    if(length(lowerBound) != n) stop("Length of lower bound of poisson random variable must be equal to n")
    ret1 <- rpois(n, lambda)
    ret <- c(mapply(max, ret1, lowerBound))
  } else {
  ret <- rpois(n, lambda)
  }
  return(ret)
}

dmypois <- nimble::nimbleRcall(
  prototype = function(
    x = double(1), #x is a matrix
    lambda = double(0),
    log = integer(0, default = 1)
  ) {},
  returnType = double(0), # outcome is a vector
  Rfun = 'mydpois'
)

rmypois <- nimble::nimbleRcall(
  prototype = function(
    n = integer(0),
    lambda = double(0),
    lowerBound = double(1, default = c(0,0)),
    includeLowerBound = integer(0, default = 0)
  ) {},
  returnType = double(1),
  Rfun = 'myrpois'
)


