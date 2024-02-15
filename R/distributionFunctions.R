
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
